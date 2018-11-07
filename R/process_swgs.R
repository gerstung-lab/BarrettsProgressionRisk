#' Run QNDAseq on bam files using previously determined parameters and output the required raw/filtered data files
#' @name runQDNAseq
#' @param bamPath Location with bam files to be processed using QDNAseq 
#' @param outputPath Location to output resulting raw and fitted read files and plots
#' 
#' @author 
#' @export
runQDNAseq<-function(bamPath, outputPath) {
  require(Biobase) 
  require(QDNAseq) 
  require(tidyverse) 
  
  if (length(list.files(bamPath, 'bam')) <= 0)
    stop("No bam files found in current directory. Exiting.")

  if (!dir.exists(outputPath))
    dir.create(outputPath, recursive = T)
  
  # get/read bin annotations
  bins <- QDNAseq::getBinAnnotations(binSize = 15)
  
  saveRDS(bins, paste(outputPath, "15kbp.rds", sep='/'))

  # write bin annotations
  pData(bins) %>%
    dplyr::mutate_at(vars(start, end), funs(as.integer)) %>%
    write_tsv(paste(outputPath,"15kbp.txt",sep='/'))
  
  # process BAM files obtaining read counts within bins
  readCounts <- QDNAseq::binReadCounts(bins, path = bamPath)
  
  pData(readCounts) %>%
    rownames_to_column(var = "sample") %>%
    dplyr::select(id = name, everything()) %>%
    write_tsv(paste(outputPath,"readCountSummary.txt",sep='/'))
  
  QDNAseq::exportBins(readCounts, paste(outputPath,"readCounts.txt",sep='/'), logTransform = FALSE)
  
  saveRDS(readCounts, paste(outputPath,"readCounts.rds",sep='/'))

  # apply filters for which bins are used, including loess residuals of calibration set
  # and blacklisted regions both taken from bin annotations
  readCountsFiltered <- QDNAseq::applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
  
  # fix for issue with copy number segmentation arising from zero count bins
  # plot median read counts as a function of GC content and mappability as an isobar plot
  pdf(paste(outputPath,"isobar.pdf",sep='/'))
  isobarPlot(readCountsFiltered)
  dev.off()
  
  # estimate the correction for GC content and mappability
  readCountsCorrected <- QDNAseq::estimateCorrection(readCountsFiltered)
  
  # output raw and fitted read counts
  features <- fData(readCountsCorrected) %>%
    as.data.frame %>%
    rownames_to_column(var = "location") %>%
    transmute(location, chrom = chromosome, start = as.integer(start), end = as.integer(end))
  
  rawReadCounts <- assayData(readCountsCorrected)$counts %>%
    as.data.frame %>%
    rownames_to_column(var = "location")
  features %>%
    dplyr::left_join(rawReadCounts, by = "location") %>%
    write_tsv(paste(outputPath,"rawReadCounts.txt",sep='/'))
  
  fittedReadCounts <- assayData(readCountsCorrected)$fit %>%
    as.data.frame %>%
    rownames_to_column(var = "location") %>%
    mutate_if(is.numeric, funs(round(., digits = 3)))
  features %>%
    dplyr::left_join(fittedReadCounts, by = "location") %>%
    write_tsv(paste(outputPath,"fittedReadCounts.txt",sep='/'))
  
  # noise plot showing relationship between the observed standard deviation in the data
  # and its read depth
  pdf(paste(outputPath,"noise.pdf",sep='/'))
  noisePlot(readCountsCorrected)
  dev.off()
}



#' @name loadSampleInformation
#' @param samples Either a filename or dataframe with appropriate information
#' @return SampleInformation object (annotated tibble)
#'
#' @author skillcoyne
#' @export
loadSampleInformation<-function(samples, path=c('NDBE','ID','LGD','HGD','IMC','OAC')) {
  if (is.character(samples)) {
    sample.info = .readFile(samples)
  } else if (is.data.frame(samples)) {
    sample.info = samples
  } else {
    stop("Filename or data frame required.")
  }

  if (length(which(c('Sample','Endoscopy') %in% colnames(sample.info))) < 2)
    stop("Sample information requires at least two columns: 'Sample' and 'Endoscopy'. Sample should be unique text identifying the sample (matching the samples in your data files), and 'Endoscopy' should be a date or integer value indicating the (descending) order in which the endoscopy was performed. Multiple samples may belong to a single endoscopy.")

  colnames(sample.info) = .titleCase(colnames(sample.info))
  
  exp_cols = c('Pathology','GEJ.Distance', 'P53 IHC')
  cols_found = sapply( exp_cols, function(x) x %in% colnames(sample.info) )
    
  if (length(which(!cols_found)) > 0)
    warning(paste0('Missing expected columns from sample information: ',paste(names(which(!cols_found)), collapse=', '), ". Recommendations will be based on predicted risks and: ",paste(names(which(cols_found)), collapse=', ')) )
  
  sample.info$Endoscopy = factor(sample.info$Endoscopy, ordered=T)
  
  pathCol = grep('Pathology',colnames(sample.info))
  if (length(pathCol) > 0) 
    sample.info[[pathCol]] = factor(sample.info[[pathCol]], levels=path, labels=path, ordered=T)
  
  p53Col = grep('P53',colnames(sample.info))
  if (length(p53Col) > 0) 
    sample.info[[p53Col]] = factor(sample.info[[p53Col]], levels=c(0,1), labels=c('Normal','Aberrant'), ordered=T)

  class(sample.info) <- c('SampleInformation', class(sample.info))
  return(sample.info)
}



#' Subtract arm-level values from smaller tile sizes. Used only if tileSegments has been run more than once with 'arms' as one tile size.
#' @name subtractArms
#' @param segments Matrix from tileSegments
#' @param arms Matrix from tileSegments (generally using 'arms' size)
#' @return matrix of values 
#'
#' @author skillcoyne
#' @export
subtractArms<-function(segments, arms) {

  if (is.null(segments) | is.null(arms))
    stop("Two matrices required")

  if (nrow(segments) != nrow(arms))
    stop(paste("Segment matrix cannot be adjusted by an arm matrix with different numbers of samples:", nrow(segments), ",", nrow(arms)))

  seg.loc = get.loc(segments)
  arm.loc = get.loc(arms)

  armsDF = makeGRangesFromDataFrame(arm.loc)
  segDF = makeGRangesFromDataFrame(seg.loc)

  tmp = segments
  # subtract arms from 5e6 and merge both (per sample)
  ov = findOverlaps(armsDF, segDF)
  for (hit in unique(queryHits(ov))) {
    cols = subjectHits(ov)[which(queryHits(ov) == hit)]
    for (i in 1:nrow(tmp)) {
      tmp[i,cols] = tmp[i,cols] - arms[i,hit]
    }
  }
  mergedDf = cbind(tmp, arms)
  return(mergedDf)
}

#' This function runs the copynumber::pcf algorithm to segment the sWGS data.
#' renamed from 'binSWGS'
#' @name segmentRawData
#' @param raw.data Data frame of raw read counts (file name is also valid)
#' @param fit.data Data frame of fitted read values (file name is also valid)
#' @param blacklist qDNAseq_blacklistedRegions (defaults to file provided in package)
#' @param min.probes minimum number of probes per segment DEF=67  (~1Mb)
#' @param gamma2 gamma adjustment for pcf DEF=250
#' @param cutoff is the residual value cutoff for QC DEF=0.015
#' @param logTransform DEF=F
#' @return list of objects:
#' 'seg.vals'=segmented samples that have passed QC, 'residuals'=data frame of per-sample residuals, 'prepped.data'=adjusted raw values, 'seg.plots'=list of per-sample genome-wide plots, 'genome.coverage'=calculated genome coverage, 'failedQC'=segmented samples that have failed QC
#'
#' @author skillcoyne
#' @export
segmentRawData<-function(info, raw.data, fit.data, blacklist=NULL, min.probes=67, gamma2=250, cutoff=0.015, logTransform=F, cache.dir=getcachedir(), build='hg19', verbose=T) {
  if (!'SampleInformation' %in% class(info))
    stop("SampleInformation object from loadSampleInformation(...) required")

  if (is.null(blacklist) | !is.data.frame(blacklist)) 
    blacklist = readr::read_tsv(system.file("extdata", "qDNAseq_blacklistedRegions.txt", package="BarrettsProgressionRisk"), col_names=T, col_types=cols(col_character(), col_integer(), col_integer()))

  if (is.character(raw.data) & is.character(fit.data)) {
    raw.data = readr::read_tsv(raw.data, col_names=T, col_types = cols('chrom'=col_character()))
    fit.data = readr::read_tsv(fit.data, col_names=T, col_type = cols('chrom'=col_character()))
  } else if (is.data.frame(fit.data) & !is.tibble(fit.data)) {
    raw.data = as_tibble(raw.data)
    fit.data = as_tibble(fit.data)
  } else {
    stop("raw and fit data must be provided as data frames with columns: location, chrom, start, end followed by sample column(s).")
  }

  chr.info = chrInfo(build=build)
  
  chrCol = grep('chr',colnames(fit.data))
  startCol = grep('start',colnames(fit.data))
  fit.data[[chrCol]] = factor(fit.data[[chrCol]], levels=levels(chr.info$chr), ordered=T)
  raw.data[[chrCol]] = factor(raw.data[[chrCol]], levels=levels(chr.info$chr), ordered=T)
  
  fit.data = fit.data %>% dplyr::arrange(fit.data[[chrCol]], fit.data[[startCol]])
  raw.data = raw.data %>% dplyr::arrange(raw.data[[chrCol]], raw.data[[startCol]])
    
  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)

  smps = intersect(info$Sample, colnames(fit.data)[countCols])
  if (length(smps) != length(info$Sample)) {
    warning(paste0("SampleInformation object and fit/raw data do not have the same samples. Using only samples from SampleInformation (n=", length(smps), ")."))
    fit.data = fit.data[,c(grep('loc|feat|chr|start|end', colnames(fit.data),value=T), smps) ]
    raw.data = raw.data[,c(grep('loc|feat|chr|start|end', colnames(raw.data),value=T), smps) ]
  }

  prepped = .prepRawSWGS(raw.data,fit.data,blacklist,logTransform)
  data = prepped$data

  sdevs = prepped$sdevs
  sdev = exp(mean(log(sdevs[!is.na(sdevs)])))
  if(verbose) message(paste('sdev=',signif(sdev,3),sep=''))

  if (ncol(data) < 4) { # Single sample
      if (verbose) message(paste("Segmenting single sample gamma=",round(gamma2*sdev,2)))
      res = copynumber::pcf( data=data, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F )
      colnames(res)[grep('mean', colnames(res))] = colnames(raw.data)[countCols]
      res$sampleID = NULL
  } else { # for most we have multiple samples
      message(paste("Segmenting", (ncol(data)-2), "samples gamma=",round(gamma2*sdev,2)))
      res = copynumber::multipcf( data=data, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F )
  }
  tmp.seg = tempfile("segments.",cache.dir,".Rdata")
  save.image(file=tmp.seg)

  if (!is.null(min.probes) & !is.na(min.probes)) {
    probes = which(res$n.probes < min.probes)
    if (verbose) message(paste(round(length(probes)/nrow(res), 3), 'segments with fewer than', min.probes, '"probes"'))
  }
  if (length(probes) > 0) res = res[-probes,]

  resids = .calculateSegmentResiduals(res, data, verbose=verbose)
  resids = resids[which(!is.na(sdevs))]
  
  coverage = round(sum(with(res, end.pos-start.pos))/chr.info[22,'genome.length'],3)
  if (verbose) message(paste(coverage, 'of the genome covered by segments.'))
  
  # --- Plot segmented data
  plist = list()
  fit.data = prepped$fit.data
  good.bins = prepped$good.bins
  window.depths.standardised = prepped$window.depths.standardised
  if (verbose) message('Plotting segmented data.')
  for(col in which(!is.na(sdevs))) {
    p = .plotSegmentedGenome(fitted = fit.data[good.bins,c(1:4,4+col)], segmented = res[,c(1:5,5+col)], window.depths.std = window.depths.standardised[good.bins,col,drop=F]) + labs(title=colnames(res)[5+col])

    med = median(res[,5+col])
    std = sd(res[,5+col])*2
    p = p + geom_hline(yintercept = c(med-std,med+std), color='grey')

    plist[[colnames(res)[5+col]]] = p
  }

  # Get mean(var(MAD(segments))) per sample
  pvr = .per.sample.residual.variance(resids)
  pvr$Pass = pvr$varMAD_median <= cutoff
  
  qcsamples = as.character(subset(pvr, Pass)$sample)
  if (verbose) message(paste(length(qcsamples), '/', nrow(pvr), ' samples passed QC.', sep=''))

  passedQC = res[,c(1:5,grep(paste(qcsamples,collapse='|'), colnames(res)))]
  failedQC = res[,unique(c(1:5,grep(paste(qcsamples,collapse='|'), colnames(res), invert=T)) )]
  if (ncol(failedQC) <= 5) failedQC = NULL

  # There's a better way to do objects/classes in R, need to spend some time with it
  swgsObj = list('seg.vals'=passedQC, 'residuals'=pvr, 'segment.residual.MSE'=resids, 'prepped.data'=data, 'seg.plots'=plist, 
             'genome.coverage'=coverage, 'failedQC'=failedQC, 'temp.file'=tmp.seg, 'cv.plot'=prepped$cv.plot, 'chr.build.info'=build, 'sample.info'=info)
  class(swgsObj) <- c('SegmentedSWGS', class(swgsObj))
  save(swgsObj, file=tmp.seg)

  return(swgsObj)
}


#' Per-sample complexity score.
#' @name scoreCX
#' @param df data frame of genome wide values
#' @param MARGIN columns or rows (see 'apply')
#' @return cx value
#'
#' @author skillcoyne
#' @export
scoreCX <- function(df, MARGIN) {
  cx = apply(df, MARGIN, function(x) {
    length(which(x >= mean(x,na.rm=T)+sd(x,na.rm=T)*2 | x <= mean(x,na.rm=T)-sd(x,na.rm=T)*2))
  })
  return(cx)
}


#' Presumes that the segmented data has already been filtered for num of probes, SD etc
#' renames 'tile.segmented.data'
#' @name tileSegments
#' @param SegmentedSWGS object from segmentRawData
#' @param size tile size to use across the genome, DEF=5e6, options are integer values or "arms"
#' @param build genome build DEF=hg19
#' @return matrix of weighted mean segmented values per tile
#'
#' @author skillcoyne
#' @export
tileSegments<-function(swgsObj, size=5e6, verbose=T) {
  require(tibble)
  
  if (class(swgsObj)[1] != 'SegmentedSWGS')
    stop("SegmentedSWGS object missing")
  
  if (!is.numeric(size) & size != 'arms')
    stop("Size must be numeric, or 'arms'")

  data = swgsObj$seg.vals
  resids = swgsObj$segment.residual.MSE
  if (!is.tibble(data)) data = as_tibble(data)
  
  #mse = mse[,intersect(colnames(data), colnames(mse))]
  
  descCols = sort(union(grep('chr|arm|start|end|probes', colnames(data), ignore.case=T), which(!sapply(data, is.numeric))))
  dataCols = c((descCols[length(descCols)]+1):ncol(data))
  
  chrCol = grep('chr',colnames(data),ignore.case=T, value=T)
  armCol = grep('arm',colnames(data),ignore.case=T,value=T)
  startPos = grep('start',colnames(data),ignore.case=T,value=T)
  endPos = grep('end',colnames(data),ignore.case=T,value=T)

  data = data[which(!data[[chrCol]] %in% c('X','Y')),]
  x1 = data[,c(chrCol, startPos, endPos, colnames(data)[dataCols])]

  tiles = .tile.genome(size, chrInfo(build=swgsObj$chr.build.info), allosomes=length(which(grepl('X|Y', unique(data[[chrCol]])))) > 0)
  gr = GenomicRanges::makeGRangesFromDataFrame(x1, keep.extra.columns=T, start.field=startPos, end.field=endPos  )
  #mseGR = GenomicRanges::makeGRangesFromDataFrame(mse[,c(chrCol,startPos,endPos,colnames(data)[dataCols])], keep.extra.columns=T, start.field=startPos, end.field=endPos  )

  mergedDf = (do.call(rbind, lapply(tiles, function(tile) {
    cbind('chr'=as.character(seqnames(tile)), as.data.frame(ranges(tile))[1:2])
  }) ))
  mergedDf[colnames(data)[dataCols]] = NA
  errorDf = mergedDf

  meanSegs = c()
  ov = GenomicRanges::findOverlaps(tiles, gr) # Each tile is a chromosome, so the overlaps are only per chromosome here
  for (chr in unique(queryHits(ov))) { # for each chromosome get overlaps
    currentChr = tiles[[chr]]
    curov = GenomicRanges::findOverlaps(currentChr, gr)

    for (i in unique(queryHits(curov))) {
      bin = currentChr[i]

      segments = gr[subjectHits(curov[queryHits(curov) == i])]

      weights = sapply(as(segments,'GRangesList'),function(r) width(pintersect(bin, r))/width(bin) )
      
      rows = with(mergedDf, which( chr==as.character(seqnames(bin)) & start == start(bin) & end == end(bin)))
      if (length(segments) > 1 & verbose)
        message(paste("chr", chr, "bin", bin, "has", length(segments), "matching segments"))

      meanSegs = c(meanSegs, length(segments))

      # weight means by the coverage of the bin
      values = apply(GenomicRanges::elementMetadata(segments), 2, weighted.mean, w=weights, na.rm=T)

      vMSE = apply(t(do.call(rbind,lapply(resids, function(sample) {
        sapply(sample[subjectHits(curov[queryHits(curov) == i])], sd)
      }))),2,weighted.mean,weights)

      #      vMSE = apply(GenomicRanges::elementMetadata(segMSE), 2, weighted.mean, w=weights, na.rm=T)
      mergedDf[rows, names(values)] = values
      errorDf[rows, names(vMSE)] = vMSE
    }
  }

  message(paste("Mean number of CN segments per genome bin:", round(mean(meanSegs, na.rm=T), 2), "median:", round(median(meanSegs, na.rm=T), 2)))

  rownames(mergedDf) = paste(mergedDf$chr, ':', mergedDf$start, '-', mergedDf$end, sep='')
  rownames(errorDf) = paste(errorDf$chr, ':', errorDf$start, '-', errorDf$end, sep='')
  mergedDf[,c(1:3)] = NULL
  errorDf[,c(1:3)] = NULL

  # deal with NA values
  mergedDf = apply(mergedDf,2, function(x) {
    x[is.na(x)] = median(x,na.rm=T)
    return(x)
  })

  errorDf = apply(errorDf,2, function(x) {
    x[is.na(x)] = median(x,na.rm=T)
    return(x)
  })
  
  return(list('tiles'=t(mergedDf), 'error'=t(errorDf)))
}


# Not used, but available if wanted now.
.logTransform<-function(x,inf=F) {
  if (is.matrix(x) | is.data.frame(x)) {
    x[] = apply(x, 2, logTV)
  } else {
    x = logTV(x)
  }
  return(x)
}

# Set up the tiles for a genome
.tile.genome<-function(tile.w=5e6, chr.info=chrInfo(build='hg19'), allosomes=F) {
  require(GenomicRanges)

  chr.info$chrom = sub('^chr', '', chr.info$chrom)

  chrs = c(1:22)
  if (allosomes) chrs = c(1:22, 'X','Y')
  chr.info = subset(chr.info, chrom %in% chrs)
  chr.info$chrom = factor(chr.info$chrom, levels=chrs, ordered=T)

  chr.info$start = 1

  genome = GenomicRanges::makeGRangesFromDataFrame(chr.info, seqnames.field = 'chrom', end.field='chr.length')
  if (is.numeric(tile.w)) {
    tiles = tile(genome, width=tile.w)
  } else if (tile.w == 'arms') {
    parms = as.data.frame(chr.info %>% rowwise %>% dplyr::summarise(
      'chr'=chrom,
      'start'=1, 'end'=chr.cent-cent.gap,
      'arm'='p'))
    qarms = as.data.frame(chr.info %>% rowwise %>% dplyr::summarise(
      'chr'=chrom,
      'start'=chr.cent+cent.gap, 'end'=chr.length,
      'arm'='q'))
    genome = makeGRangesFromDataFrame(rbind(parms,qarms), seqnames.field = 'chr', keep.extra.columns = T)
    tiles = lapply( levels(seqnames(genome)), function(seq) genome[seqnames(genome)==seq] )
    tiles = GRangesList(tiles)
    if (length(tiles) != length(unique(seqnames(genome)))) stop("Tile and genome length doesn't match")
  } else if (tile.w == 'chr') {
    tiles = lapply(levels(seqnames(genome)), function(seq) genome[seqnames(genome)==seq] )
    tiles = GRangesList(tiles)
    if (length(tiles) != length(unique(seqnames(genome)))) stop("Tile and genome length doesn't match")
  }
  return(tiles)
}

# Calculate the residuals per segment between the fitted values from QDNAseq, and the segments post-pcf
# DISTANCE between segment and points, then variance across the distances...
.calculateSegmentResiduals<-function(calcSegments, observedCN, verbose=T) {
  if (verbose) message("Calculating variance of residuals per segment")
  cols = ncol(observedCN)-2
  resids = list()
  for (i in 1:cols) {
    sample.name = colnames(observedCN)[-c(1:2)][i]
    pred.seg = calcSegments[,c(1:5,(5+i))]

    resvar = lapply(1:nrow(pred.seg), function(j) {
      seg = pred.seg[j,]
      rows = which(observedCN$start >= seg[['start.pos']] & observedCN$start <= seg[['end.pos']] )
      (seg[[sample.name]]-observedCN[rows,sample.name])
    })
    resids[[sample.name]] = resvar
  }
  return(resids)
}


# Calculate per-sample residual variance from the list generated in .calculateSegmentResiduals
.per.sample.residual.variance<-function(segment.residuals) {
  if (is.null(segment.residuals) | !is.list(segment.residuals)) {
    stop("List of segment residuals required")
  }

  cols = c('samplename','varMAD','n.segs')
  res.variance = (data.frame(matrix(ncol=length(cols), nrow=0, dimnames=list(c(),cols))))

  var.resids = lapply(segment.residuals, function(sample) {
    do.call(rbind.data.frame, lapply(sample, function(y) {
      cbind('varMAD'=var(y[y<mad(y) & y>-mad(y)]))
    }))
  })
  for (sample in names(var.resids)) {
    n.segs = length(segment.residuals[[sample]])

    msd = var.resids[[sample]] %>% dplyr::summarise_all(funs(median,sd) )
    q1 = var.resids[[sample]] %>% dplyr::summarise_all(funs(Q1=quantile),probs=0.25 )
    q3 = var.resids[[sample]] %>% dplyr::summarise_all(funs(Q3=quantile),probs=0.75 )

    res.variance = rbind(res.variance, cbind(sample,msd,q1,q3,n.segs))
  }
  colnames(res.variance)[2:5] = paste('varMAD_',colnames(res.variance)[2:5],sep='')

  return(res.variance)
}

# preps data for segmentation
.prepRawSWGS<-function(raw.data,fit.data,blacklist, logTransform=F) {
  if (ncol(blacklist) < 3)
    stop('Blacklisted regions missing or incorrectly formatted.\nExpected columnes: chromosome start end')

  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)
  infoCols = grep('loc|feat|chr|start|end', colnames(fit.data))

  rows = intersect(fit.data[[1]], raw.data[[1]])
  #rows = which(fit.data[[1]] %in% raw.data[[1]])
  fit.data = fit.data[which(fit.data[[1]] %in% rows),]
  raw.data = raw.data[which(raw.data[[1]] %in% rows),]

  window.depths = as.vector(as.matrix(raw.data[,countCols]))/as.vector(as.matrix(fit.data[,countCols]))

  # QDNAseq does this in 'correctBins' but we don't use that method so added here
  negs = which(as.vector(as.matrix(fit.data[,countCols])) <= 0)
  if (length(negs) > 0) window.depths[negs] = 0

  window.depths = matrix(window.depths, ncol=length(countCols))

  chr.info = chrInfo()
  plist = list()
  
  df = cbind(fit.data[,infoCols], window.depths)
  colnames(df)[-infoCols] = colnames(fit.data)[countCols]
  df = base::merge(df, chr.info[,c('chr','chr.length')], by.x='chrom',by.y='chr',all.x=T)

  df = reshape::melt(df, measure.vars=colnames(fit.data)[countCols])
  pp = ggplot(df, aes(x=1:chr.length)) + ylim(c(0, quantile(df$value, probs=0.75, na.rm=T)*2)) +
    facet_grid(variable~chrom, space='free', scales='free') +
    geom_point( aes(start, value), color='darkred', alpha=.4) +  
    labs(title='Coverage', x='Chromosomes', y="corrected depth/15KB") + 
    theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0,'lines'))

  message(paste(nrow(blacklist), "genomic regions in the exclusion list."))

  fit.data$in.blacklist = F
  for(r in 1:nrow(blacklist)) {
    fit.data$in.blacklist[ fit.data$chrom == blacklist$chromosome[r] &
                             fit.data$start >= blacklist$start[r] &
                             fit.data$end <= blacklist$end[r] ] = T
  }
  message(paste("# blacklisted probes = ",sum(fit.data$in.blacklist), ' (',round(sum(fit.data$in.blacklist)/nrow(fit.data),2)*100,'%)',sep=""))

  if (sum(fit.data$in.blacklist) <= 0)
    warning("No probes excluded from the blacklist.")

  if (length(countCols) == 1)
    window.depths = as.data.frame(window.depths)

  window.depths.standardised = as.data.frame(window.depths[which(!fit.data$in.blacklist),])
  if (logTransform)
    window.depths.standardised = log2( window.depths.standardised+abs(min(window.depths.standardised, na.rm=T))+1 )

  fit.data = fit.data[!fit.data$in.blacklist,-ncol(fit.data)]
  sdevs = sapply(c(1:length(countCols)), function(s) {
    getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25 )
  })
  sdevs[sdevs==0] = NA

  good.bins = which(!is.na(rowSums(as.data.frame(window.depths.standardised[,!is.na(sdevs)]))))

  data = cbind(fit.data[good.bins,c('chrom','start')],window.depths.standardised[good.bins,!is.na(sdevs)])
  colnames(data)[-c(1:2)] = colnames(raw.data)[countCols]

  return(list('data'=data,'sdevs'=sdevs, 'good.bins'=good.bins, 'window.depths.standardised'=window.depths.standardised, 'fit.data'=fit.data, 'cv.plot'=pp))
}
