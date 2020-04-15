#' Run QNDAseq on bam files using previously determined parameters and output the required raw/filtered data files
#' @name runQDNAseq
#' @param bamPath Location with bam files to be processed using QDNAseq 
#' @param outputPath Location to output resulting raw and fitted read files and plots
#' 
#' @author skillcoyne
#' @export
runQDNAseq<-function(bam=NULL,path=NULL,outputPath=NULL, minMapQ=37, binsize=50) {
  require(Biobase) 
  require(QDNAseq) 
  require(tidyverse) 
  
  if ((!is.null(bam) && !file.exists(bam)) || rev(unlist(strsplit(basename(bam), '\\.')))[1] != 'bam' ) {
    error = paste0('BAM file "', bam, '" does not exist.')
  } else if (is.null(bam) && length(list.files(path, 'bam')) <= 0) {
    error = paste(error, "No bam files found in current directory. Exiting.", sep='\n')
  }

  if (!is.null(error)) 
    stop(error)

  if (binsize != 50) 
    warning("Internal model was generated using data processed with a 50kb bin size for QDNAseq. Using a different bin size is not recommended.")
  
  if (!dir.exists(outputPath))
    dir.create(outputPath, recursive = T)
  
  # get/read bin annotations
  bins <- QDNAseq::getBinAnnotations(binSize = binsize)
  
  saveRDS(bins, paste(outputPath, paste0(binsize,"kbp.rds"), sep='/'))

  # write bin annotations
  pData(bins) %>%
    dplyr::mutate_at(vars(start, end), list(as.integer)) %>%
    write_tsv(paste(outputPath,paste0(binsize,"kbp.txt"),sep='/'))
  
  # process BAM files obtaining read counts within bins
  if (!is.null(bam)) {
    message(paste0('Binning read counts (min mapQ=',minMapQ,') from ', bam))
    readCounts <- QDNAseq::binReadCounts(bins, bamfiles=bam, minMapq=minMapQ)
  } else {
    message(paste0('Binning read counts from BAM files in ', path))
    readCounts <- QDNAseq::binReadCounts(bins, path=path, minMapq=minMapQ)
  }
  
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
    mutate_if(is.numeric, list(~round(., digits = 3)))
  features %>%
    dplyr::left_join(fittedReadCounts, by = "location") %>%
    write_tsv(paste(outputPath,"fittedReadCounts.txt",sep='/'))
  
  # noise plot showing relationship between the observed standard deviation in the data
  # and its read depth
  pdf(paste(outputPath,"noise.pdf",sep='/'))
  noisePlot(readCountsCorrected)
  dev.off()
}

#' Set up the sample information required for analysis.
#' @name loadSampleInformation
#' @param samples Either a filename or dataframe with appropriate information
#' @param path Ordered list of pathology abbreviations (OPT)
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
  
  col_regex = 'Pathology|GEJ(\\.| )Distance|P53(\\.| )IHC'
  cols_found = grep(col_regex, colnames(sample.info), value=T, ignore.case = T)

  if (length(cols_found) < 3) {
    message = paste0("Missing expected columns from sample information. Recommendations will be based on predicted risks")
    if (length(cols_found) > 0) message = paste0(message, ' and: ', paste(cols_found,collapse=','))
    warning(message)
  }

  endo.date<-function(endo) {
    if (is.character(endo)) {
      strings = unlist(strsplit(endo, '-|/'))
      sep = ifelse(grepl("-", endo), "-", "/")
      dates = ifelse(as.integer(strings[1]) > 1900, paste(c("%Y", "%m", "%d"), collapse = sep), paste(c("%d", "%m", "%Y"), collapse = sep))
      parse_date(endo, format=dates)
    } else if (is.numeric(endo)) {
      #def = as.Date('2001/01/01')
      return(endo)
    } else { 
      as.Date(endo)
    }
  }
  
  # Don't change the factor if it's already done.
  gej_col = grep('GEJ(\\.| )Distance', colnames(sample.info), value=T)
  if (length(gej_col) > 0 && !is.factor(sample.info[[gej_col]])) {
    gej.dist<-function(gej) {
        funct = case_when(
          length(which(is.na(as.numeric(gej)))) <= 0 ~ 'as.numeric',
          length(which(is.na(as.character(gej)))) <= 0 ~ 'as.character'
        )
        sapply(gej, eval(funct))
      }
  
    sample.info = sample.info %>% 
      dplyr::mutate_at(dplyr::vars(!!gej_col), list(~gej.dist(.)))  
    
    levels = sort(sample.info %>% dplyr::select(!!gej_col) %>% distinct %>% pull)
  
    sample.info = sample.info %>% 
      dplyr::mutate_at(dplyr::vars(!!gej_col), list(~factor(., levels=levels, ordered=T)))
  }  
  
  p53_col = grep('P53', colnames(sample.info), value=T)
  if (length(p53_col) > 0 && 
      length(grep('Normal|Aberrant', dplyr::select(sample.info, !!p53_col) %>% distinct %>% pull)) > 0) {
    sample.info = sample.info %>% 
      dplyr::mutate_at(dplyr::vars(matches('P53')), list(~factor(., levels=c('Normal','Aberrant'), ordered=T)))
  } else {
    sample.info = sample.info %>% dplyr::mutate_at(dplyr::vars(dplyr::matches('P53')), list(~factor(., levels=c(0,1), labels=c('Normal','Aberrant'), ordered=T)))
  }
  
  sample.info = sample.info %>% rowwise %>% dplyr::mutate(Endoscopy = endo.date(Endoscopy))

  sample.info = sample.info %>% 
    dplyr::mutate_at(dplyr::vars(dplyr::matches('Pathology')), list(~factor(., levels=path, ordered=T))) 

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


qdna.to.probes<-function(kb) {
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  if (!is.wholenumber(kb)) stop('Missing integer value for kb bin size used in QDNAseq segmentation')
  round(1e6/(kb*1000))
}

#' This function runs the copynumber::pcf algorithm to segment the sWGS data.
#' renamed from 'binSWGS'
#' @name segmentRawData
#' @param raw.data Data frame of raw read counts (file name is also valid)
#' @param fit.data Data frame of fitted read values (file name is also valid)
#' @param blacklist qDNAseq_blacklistedRegions (defaults to file provided in package)
#' @param kb QDNAseq bin size
#' @param gamma2 gamma adjustment for pcf DEF=250
#' @param cutoff is the residual value cutoff for QC DEF=0.015
#' @param norm normalize pcf, DEF=T
#' @param logTransform DEF=F
#' @return list of objects:
#' 'seg.vals'=segmented samples that have passed QC, 'residuals'=data frame of per-sample residuals, 'prepped.data'=adjusted raw values, 'seg.plots'=list of per-sample genome-wide plots, 'genome.coverage'=calculated genome coverage, 'failedQC'=segmented samples that have failed QC
#'
#' @author skillcoyne
#' @export
segmentRawData<-function(info, raw.data, fit.data, blacklist=readr::read_tsv(system.file("extdata", "qDNAseq_blacklistedRegions.txt", package="BarrettsProgressionRisk"), col_names=T, col_types='cii'), gamma2=250, kb=50, cutoff=0.008, multipcf=T, logTransform=F, cache.dir=getcachedir(), build='hg19', verbose=T) {
  #intPloidy=F  # This wasn't terribly useful. Leaving the code in place for the moment but setting it to false by default.
  #if (intPloidy & cutoff < 0.03) cutoff = cutoff*2

  if (cutoff %% 1 == 0) cutoff = cutoff + 0.0001

  min.probes = qdna.to.probes(kb)

  if (!'SampleInformation' %in% class(info))
    stop("SampleInformation object from loadSampleInformation(...) required")

  if (is.character(raw.data) & is.character(fit.data)) {
    raw.data = readr::read_tsv(raw.data, col_names=T, col_types = cols('chrom'=col_character()))
    fit.data = readr::read_tsv(fit.data, col_names=T, col_type = cols('chrom'=col_character()))
  } else if (is.data.frame(fit.data)) {
    raw.data = as_tibble(raw.data)
    fit.data = as_tibble(fit.data)
  } else {
    stop("raw and fit data must be provided as data frames with columns: location, chrom, start, end followed by sample column(s).")
  }

  qkb = raw.data %>% dplyr::mutate(qkb = (end-start)/1e3) %>% dplyr::summarise(qkb=round(mean(qkb))) %>% pull
  
  if (qkb != BarrettsProgressionRisk:::be_model$qdnaseq_kb)
    warning(paste0("Model requires QDNAseq bin size of ",BarrettsProgressionRisk:::be_model$qdnaseq_kb,'kb, data processed at ',qkb,'kb, predictions will be inaccurate.'))
  
  if (qkb != kb)
    stop(paste0("QDNAseq data was processed at bin size ", qkb, "kb, kb parameter is set to ", kb,'kb.'))
  
  chr.info = chrInfo(build=build)
  
  fit.data = fit.data %>% dplyr::mutate_at(vars(dplyr::matches('chr')),list(~factor(.,levels=levels(chr.info$chr), ordered=T))) %>%
    dplyr::arrange_at( vars(dplyr::matches('chr'), dplyr::matches('start')), list() )
  
  raw.data = raw.data %>% mutate_at(vars(dplyr::matches('chr')), list(~factor(.,levels=levels(chr.info$chr), ordered=T))) %>% 
    dplyr::arrange_at( vars(dplyr::matches('chr'), dplyr::matches('start')), list() )
  
  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)
  if (length(countCols) == 1) {
    colnames(raw.data)[countCols] = info$Sample
    colnames(fit.data)[countCols] = info$Sample
  }
  
  smps = intersect(info$Sample, colnames(fit.data)[countCols])
  if (length(smps) <= 0)
    stop(paste0('No matching samples in SampleInformation object and data files. Stopping.'))
  
  if ( (length(smps) != length(info$Sample)) | length(countCols) != length(info$Sample) ) {
    warning(paste0("SampleInformation object and fit/raw data do not have the same samples. Using only samples from SampleInformation (n=", length(smps), ")."))
    fit.data = fit.data[,c(grep('loc|feat|chr|start|end', colnames(fit.data),value=T), smps) ]
    raw.data = raw.data[,c(grep('loc|feat|chr|start|end', colnames(raw.data),value=T), smps) ]
    
    countCols = countCols[which(colnames(fit.data)[countCols] %in% smps)]
  }

  prepped = prepRawSWGS(raw.data,fit.data,blacklist,F)
  data = prepped$data

  raw.variance = data %>% summarise_at(vars(-chrom, -start), list(~sd(.,na.rm=T))) 
  # TODO load St.Dev data for all training 'prepped' data into the sysdata file
  if (length(which(raw.variance > 0.21) > 0))
      warning( paste(paste(names(raw.variance[which(raw.variance > 0.15)]), collapse=', '), 'raw values have a std.dev. > 0.21. Rerun QDNAseq with a larger binsize.') )

  sdevs = prepped$sdevs
  sdev = exp(mean(log(sdevs[!is.na(sdevs)])))
  if(verbose) message(paste('sdev=',signif(sdev,3),sep=''))

#  segs = list()
  ## TODO Using IntPloidy we may want to only pcf a single sample at a time
  if (ncol(data) < 4) { # Single sample
      if (verbose) message(paste("Segmenting single sample gamma=",round(gamma2*sdev,2)))
      res = copynumber::pcf( data=data, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F, assembly=build)
      colnames(res)[grep('mean', colnames(res))] = colnames(raw.data)[countCols]
      res$sampleID = NULL
#      segs[smps[1]] = res
  } else if (!multipcf) {
    if (verbose) message(paste("Segmenting multiple samples individually gamma=",round(gamma2*sdev,2)))
    
    segs = lapply(smps, function(s) {
      if (verbose) message(s)
      tmp = data %>% dplyr::select(chrom,start,dplyr::matches(s))
      res = copynumber::pcf( data=tmp, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F, assembly=build)
      colnames(res)[grep('mean', colnames(res))] = s
      res$sampleID = NULL
      return(res)
    })
    names(segs) = smps
    res = do.call(bind_rows, segs) %>% arrange(chrom, start.pos)
    
  } else { # for most we have multiple samples
      message(paste("Segmenting", (ncol(data)-2), "samples gamma=",round(gamma2*sdev,2)))
      res = copynumber::multipcf( data=data, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F, assembly=build)
  }

  tmp.seg = tempfile("segments.",cache.dir,".Rdata")
  save.image(file=tmp.seg)

  prb<-function(x) {
    x = na.omit(x)
    probes = which(x$n.probes < min.probes)
    round(length(probes)/nrow(x),3)
  }
  prb.ratio = sapply(smps, function(s)  prb(res %>% dplyr::select(chrom, n.probes, dplyr::matches(s))) )
  if (verbose) {
    message(paste0('Per sample ratio of segments with fewer than ', min.probes, ' "probes" - '))
    print(prb.ratio)
  }
    
  probes = which(res$n.probes < min.probes)
  if (length(probes) > 0) res = res[-probes,]

  resids = .calculateSegmentResiduals(res, data, verbose=verbose)
  resids = resids[which(!is.na(sdevs))]
  
  #if (intPloidy) res = res %>% dplyr::mutate_at(vars(info$Sample), list( ~round(.,1) ))

  cvg<-function(x) {
    x = na.omit(x)
    round(sum(as.numeric(with(x, end.pos-start.pos)),na.rm=T)/(chr.info %>% filter(chr == 22) %>% dplyr::select(genome.length) %>% pull),3)
  }
  
  coverage = sapply(smps, function(s)  cvg(res %>% dplyr::select(chrom, start.pos, end.pos, dplyr::matches(s))) )
  #coverage = round(sum(as.numeric(with(res, end.pos-start.pos)),na.rm=T)/chr.info[22,'genome.length'],3)
  if (verbose) message(paste(round(mean(coverage),2), 'of the genome covered by segments.'))
  
  # --- Plot segmented data
  plist = tryCatch({
    plist = list()
    fit.data = prepped$fit.data
    good.bins = prepped$good.bins
    window.depths.standardised = prepped$window.depths.standardised
    if (verbose) message('Plotting segmented data.')
    for(col in which(!is.na(sdevs))) {
      p = .plotSegmentedGenome(fitted = fit.data[good.bins,c(1:4,4+col)], segmented = na.omit(res[,c(1:5,5+col)]), window.depths.std = window.depths.standardised[good.bins,col,drop=F]) + labs(title=colnames(res)[5+col])
  
      med = median(res[,5+col],na.rm=T)
      std = sd(res[,5+col],na.rm=T)*2
      p = p + geom_hline(yintercept = c(med-std,med+std), color='grey')
  
      plist[[colnames(res)[5+col]]] = p
    }
    plist
  }, error = function(e) {
    warning(paste0("Failed to plot segments: ", e))
    NULL
  })

  # Get mean(var(MAD(segments))) per sample
  if (verbose) message('Calculating sample residual variance.')
  pvr = .per.sample.residual.variance(resids)
  digits = length(unlist((strsplit(unlist(strsplit(as.character(cutoff), '\\.'))[2], ''))))
  if (digits == 1) cutoff  
  pvr = pvr %>% mutate_at(vars(contains('MAD')), list(round), digits) %>% mutate(Pass = round(pvr$varMAD_median, 3) <= cutoff)
  
  qcsamples = pvr %>% filter(Pass) %>% dplyr::select(dplyr::matches('sample')) %>% pull
  if (verbose) message(paste(length(qcsamples), '/', nrow(pvr), ' samples passed QC.', sep=''))

  passedQC = res[,c(1:5,grep(paste(qcsamples,collapse='|'), colnames(res)))]
  failedQC = res[,unique(c(1:5,grep(paste(qcsamples,collapse='|'), colnames(res), invert=T)) )]
  if (ncol(failedQC) <= 5) failedQC = NULL

  # There's a better way to do objects/classes in R, need to spend some time with it
  swgsObj = list('seg.vals'=passedQC, 'residuals'=pvr, 'segment.residual.MSE'=resids, 'prepped.data'=data, 'seg.plots'=plist, 
             'genome.coverage'=coverage, 'failedQC'=failedQC, 'temp.file'=tmp.seg, 'cv.plot'=prepped$cv.plot, 'chr.build.info'=build, 'sample.info'=info, bin.size=kb)
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
  checkErr = T
  if (class(swgsObj)[1] != 'SegmentedSWGS') {
    warning("SegmentedSWGS object missing, per-tile errors and sample QC cannot be assessed.")
    checkErr = F
    data = swgsObj
    build = 'hg19'
    failed = NULL
  } else {
    data = swgsObj$seg.vals
    failed = sampleResiduals(swgsObj) %>% dplyr::filter(!Pass)
    build = swgsObj$chr.build.info
    if (nrow(failed) == nrow(sampleResiduals(swgsObj) ) )
      stop('All samples failed QC, no data available for prediction.')
  }
  
  if (!is.numeric(size) & size != 'arms')
    stop("Size must be numeric, or 'arms'")

  if (checkErr) {
    resids = swgsObj$segment.residual.MSE[as.character((swgsObj$residuals %>% dplyr::filter(Pass) %>% dplyr::select(dplyr::matches('sample')) %>% pull))] 
  }
  
  if (!is_tibble(data)) data = as_tibble(data)
  
  #mse = mse[,intersect(colnames(data), colnames(mse))]
  
  descCols = sort(union(grep('chr|arm|start|end|probes', colnames(data), ignore.case=T), which(!sapply(data, is.numeric))))
  dataCols = c((descCols[length(descCols)]+1):ncol(data))
  
  chrCol = grep('chr',colnames(data),ignore.case=T, value=T)
  armCol = grep('arm',colnames(data),ignore.case=T,value=T)
  startPos = grep('start',colnames(data),ignore.case=T,value=T)
  endPos = grep('end',colnames(data),ignore.case=T,value=T)

  # No sex chromosomes
  data = data %>% filter(!!sym(chrCol) %in% c(1:22) )
  
  tiles = .tile.genome(size, chrInfo(build=build), allosomes=length(which(grepl('X|Y', unique(data[[chrCol]])))) > 0)
  mergedDf = (do.call(rbind, lapply(tiles, function(tile) {
    cbind('chr'=as.character(seqnames(tile)), as.data.frame(ranges(tile))[1:2])
  }) ))
  mergedDf[colnames(data)[dataCols]] = NA
  errorDf = mergedDf
  
  x1 = data %>% dplyr::select(chrCol, startPos, endPos, colnames(data)[dataCols])

  gr = GenomicRanges::makeGRangesFromDataFrame(x1, keep.extra.columns=T, start.field=startPos, end.field=endPos  )
  #mseGR = GenomicRanges::makeGRangesFromDataFrame(mse[,c(chrCol,startPos,endPos,colnames(data)[dataCols])], keep.extra.columns=T, start.field=startPos, end.field=endPos  )

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
      if (length(segments) > 0 & verbose)
        message(paste("chr", chr, "bin", bin, "has", length(segments), "matching segments in", length(dataCols), 'samples.'))

      meanSegs = c(meanSegs, length(segments))

      # weight means by the coverage of the bin
      values = apply(GenomicRanges::elementMetadata(segments), 2, weighted.mean, w=weights, na.rm=T)

      if (checkErr) {
        vMSE = apply(t(do.call(rbind, lapply(resids, function(sample) {
          sapply(sample[subjectHits(curov[queryHits(curov) == i])], sd)
        }))),2,weighted.mean,weights)
      }

      #      vMSE = apply(GenomicRanges::elementMetadata(segMSE), 2, weighted.mean, w=weights, na.rm=T)
      mergedDf[rows, names(values)] = values
      if (checkErr) errorDf[rows, names(vMSE)] = vMSE
    }
  }

  if (verbose)
    message(paste("Mean number of CN segments per genome bin:", round(mean(meanSegs, na.rm=T), 2), "median:", round(median(meanSegs, na.rm=T), 2)))

  rownames(mergedDf) = paste(mergedDf$chr, ':', mergedDf$start, '-', mergedDf$end, sep='')
  mergedDf[,c(1:3)] = NULL
  if (checkErr) {
    rownames(errorDf) = paste(errorDf$chr, ':', errorDf$start, '-', errorDf$end, sep='')
    errorDf[,c(1:3)] = NULL
  }
  
  # deal with NA values
  mergedDf = apply(mergedDf,2, function(x) {
    x[is.na(x)] = median(x,na.rm=T)
    return(x)
  })

  if (checkErr) {
    errorDf = apply(errorDf,2, function(x) {
      x[is.na(x)] = median(x,na.rm=T)
      return(x)
    })
  }
  
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
    #pred.seg = na.omit(calcSegments[,c(1:5,(5+i))])

    pred.seg = na.omit(calcSegments %>% dplyr::select(chrom, arm, start.pos, end.pos, n.probes, !!sample.name))
    
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

  var.resids = lapply(segment.residuals, function(sample) {
    do.call(rbind.data.frame, lapply(sample, function(y) {
      cbind('varMAD'=var(y[y<mad(y,na.rm=T) & y>-mad(y,na.rm=T)]))
    }))
  })
  
  res.variance = tibble()
  for (sample in names(var.resids)) {
    n.segs = length(segment.residuals[[sample]])

    msd = var.resids[[sample]] %>% dplyr::summarise_all(list(~median(.,na.rm=T),~sd(.,na.rm=T)) )
    q1 = var.resids[[sample]] %>% dplyr::summarise_all(list(Q1=quantile),probs=0.25,na.rm=T )
    q3 = var.resids[[sample]] %>% dplyr::summarise_all(list(Q3=quantile),probs=0.75,na.rm=T )

    res.variance = bind_rows(res.variance, bind_cols(bind_cols(msd,q1), q3) %>% add_column(samplename = sample, .before=1))
  }
  colnames(res.variance)[2:5] = paste('varMAD_',colnames(res.variance)[2:5],sep='')

  return(res.variance)
}

# preps data for segmentation
prepRawSWGS<-function(raw.data,fit.data,blacklist = readr::read_tsv(system.file("extdata", "qDNAseq_blacklistedRegions.txt", package="BarrettsProgressionRisk"), col_names=T, col_types='cii'), plot=T,verbose=F) {

#  intPloidy=F
  
  if (ncol(blacklist) < 3)
    stop('Blacklisted regions missing or incorrectly formatted.\nExpected columnes: chromosome start end')

  rows = intersect(fit.data[[1]], raw.data[[1]])
  fit.data = fit.data[which(fit.data[[1]] %in% rows),]
  raw.data = raw.data[which(raw.data[[1]] %in% rows),]

  sortedCountCols = intersect(colnames(raw.data), colnames(fit.data))

  raw.data = raw.data %>% dplyr::select(dplyr::matches('loc|feat|chr|start|end'), sortedCountCols)
  fit.data = fit.data %>% dplyr::select(dplyr::matches('loc|feat|chr|start|end'), sortedCountCols)
  
  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)
  infoCols = grep('loc|feat|chr|start|end', colnames(fit.data))

  # raw counts adjusted by the fitted
  window.depths = as.vector(as.matrix(raw.data[,countCols]))/as.vector(as.matrix(fit.data[,countCols]))
  
  # QDNAseq does this in 'correctBins' but we don't use that method so added here
  negs = which(as.vector(as.matrix(fit.data[,countCols])) <= 0)
  if (length(negs) > 0) window.depths[negs] = 0

  window.depths = matrix(window.depths, ncol=length(countCols))

  plotlist = list()
  if (plot) {
    chr.info = chrInfo()
    
    df = cbind(fit.data[,infoCols], window.depths)
    colnames(df)[-infoCols] = colnames(fit.data)[countCols]
    df = base::merge(df, chr.info[,c('chr','chr.length')], by.x='chrom',by.y='chr',all.x=T)
    for (col in countCols) {
      tmp = reshape::melt(df[,c(infoCols, col)], measure.vars=colnames(fit.data)[col])
      p = ggplot(tmp, aes(x=1:chr.length)) + ylim(c(0, quantile(tmp$value, probs=0.75, na.rm=T)*2)) +
        facet_grid(~chrom, space='free', scales='free') +
        geom_point( aes(start, value), color='darkred', alpha=.4) +  
        labs(title=colnames(fit.data[col]), x='Chromosomes', y="corrected depth") + 
        theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0,'lines'))
      plotlist[[colnames(fit.data)[col]]] = p
    }
  }
  
  if (verbose) message(paste(nrow(blacklist), "genomic regions in the exclusion list."))

  fit.data$in.blacklist = F
  for(r in 1:nrow(blacklist)) {
    fit.data$in.blacklist[ fit.data$chrom == blacklist$chromosome[r] &
                             fit.data$start >= blacklist$start[r] &
                             fit.data$end <= blacklist$end[r] ] = T
  }
  if (verbose) message(paste("# blacklisted probes = ",sum(fit.data$in.blacklist), ' (',round(sum(fit.data$in.blacklist)/nrow(fit.data),2)*100,'%)',sep=""))

  if (sum(fit.data$in.blacklist) <= 0) warning("No probes excluded from the blacklist.")

  if (length(countCols) == 1) window.depths = as.data.frame(window.depths)

  window.depths.standardised = as.data.frame(window.depths[which(!fit.data$in.blacklist),])
  
  # transformations?  Not sure about the log transformation
  #window.depths.standardised = sqrt( window.depths.standardised+abs(min(window.depths.standardised, na.rm=T)) * 0.375 ) 
  #if (logTransform) window.depths.standardised = log2( window.depths.standardised+abs(min(window.depths.standardised, na.rm=T))+.Machine$double.xmin ) # this was the wrong way anyhow

  fit.data = fit.data[!fit.data$in.blacklist,-ncol(fit.data)]
  sdevs = sapply(c(1:length(countCols)), function(s) {
    getMad( window.depths.standardised[!is.na(window.depths.standardised[,s]),s], k=25 )
  })
  sdevs[sdevs==0] = NA

  good.bins = which(!is.na(rowSums(as.data.frame(window.depths.standardised[,!is.na(sdevs)]))))

  data = cbind(fit.data[good.bins,c('chrom','start')], window.depths.standardised[good.bins,!is.na(sdevs)])
  colnames(data)[-c(1:2)] = colnames(raw.data)[countCols[!is.na(sdevs)]]

  return(list('data'=data,'sdevs'=sdevs, 'good.bins'=good.bins, 'window.depths.standardised'=window.depths.standardised, 'fit.data'=fit.data, 'cv.plot'=plotlist))
}
