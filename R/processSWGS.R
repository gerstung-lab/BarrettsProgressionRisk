

#' Subtract arm-level values from smaller tile sizes. Used only if tileSegments has been run more than once with 'arms' as one tile size.
#' @name subtractArms
#' @param segments Matrix from tileSegments
#' @param arms Matrix from tileSegments (generally using 'arms' size)
#' @return matrix of values 
#'
#' @author skillcoyne
#' @export
subtractArms<-function(segments, arms) {
  get.loc<-function(df) {
    locs = do.call(rbind.data.frame, lapply(colnames(df), function(x) unlist(strsplit( x, ':|-'))))
    colnames(locs) = c('chr','start','end')
    locs[c('start','end')] = lapply(locs[c('start','end')], function(x) as.numeric(as.character(x)))
    locs$chr = factor(locs$chr, levels=c(1:22), ordered=T)
    locs
  }

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
  #head(mergedDf)
  #summary(mergedDf[,1])

  #plot(apply(segs,2,mean))
  #plot(apply(mergedDf, 2, mean))
  return(mergedDf)
}

#' This function runs the copynumber::pcf algorithm to segment the sWGS data.
#' renamed from 'binSWGS'
#' @name segmentRawData
#' @param raw.data Data frame of raw read counts
#' @param fit.data Data frame of fitted read values 
#' @param blacklist DEF qDNAseq_blacklistedRegions
#' @param min.probes minimum number of probes per segment DEF=67  (~1Mb)
#' @param gamma2 gamma adjustment for pcf DEF=250
#' @param cutoff is the residual value cutoff for QC DEF=0.015
#' @param logTransform DEF=F
#' @return list of objects:
#' 'seg.vals'=segmented samples that have passed QC, 'residuals'=data frame of per-sample residuals, 'prepped.data'=adjusted raw values, 'seg.plots'=list of per-sample genome-wide plots, 'genome.coverage'=calculated genome coverage, 'failedQC'=segmented samples that have failed QC
#'
#' @author skillcoyne
#' @export
segmentRawData<-function(raw.data, fit.data, blacklist=read.table(system.file("extdata", "qDNAseq_blacklistedRegions.txt", package="BarrettsProgressionRisk"), header=T, sep='\t'), min.probes=67, gamma2=250, cutoff=0.015, logTransform=F, verbose=T) {
  require(copynumber)

  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)

  prepped = .prepRawSWGS(raw.data,fit.data,blacklist,logTransform)
  data = prepped$data

  sdevs = prepped$sdevs
  sdev = exp(mean(log(sdevs[!is.na(sdevs)])))
  if(verbose) message(paste('sdev=',signif(sdev,3),sep=''))

  if (ncol(data) < 4) { # Single sample
      if (verbose) message(paste("Segmenting single sample gamma=",round(gamma2*sdev,2)))
      res = pcf( data=data, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F )
      colnames(res)[grep('mean', colnames(res))] = colnames(raw.data)[countCols]
      res$sampleID = NULL
  } else { # for most we have multiple samples
      message(paste("Segmenting", (ncol(data)-2), "samples gamma=",gamma2*sdev))
      res = multipcf( data=data, gamma=gamma2*sdev, fast=F, verbose=verbose, return.est=F )
  }

  if (!is.null(min.probes) & !is.na(min.probes)) {
    probes = which(res$n.probes < min.probes)
    if (verbose) message(paste(round(length(probes)/nrow(res), 3), 'segments with fewer than', min.probes, '"probes"'))
  }
  if (length(probes) > 0) res = res[-probes,]

  resids = .calculateSegmentResiduals(res, data, verbose=verbose)
  resids = resids[which(!is.na(sdevs))]

  chr.info = chrInfo(build='hg19')
  coverage = round(sum(with(res, end.pos-start.pos))/chr.info[22,'genome.length'],3)

  # --- Plot segmented data
  plist = list()
  fit.data = prepped$fit.data
  good.bins = prepped$good.bins
  window.depths.standardised = prepped$window.depths.standardised
  for(col in which(!is.na(sdevs))) {
    p = .plotSegmentedGenome(fitted = fit.data[good.bins,c(1:4,4+col)], segmented = res[,c(1:5,5+col)], window.depths.std = window.depths.standardised[good.bins,col,drop=F]) + labs(title=colnames(res)[5+col])

    med = median(res[,5+col])
    std = sd(res[,5+col])*2
    p = p + geom_hline(yintercept = c(med-std,med+std), color='grey')

    plist[[col]] = p
  }

  pvr = .per.sample.residual.variance(resids)
  pvr$Pass = pvr$varMAD_median <= cutoff
  
  qcsamples = as.character(subset(pvr, Pass)$sample)
  if (verbose) message(paste(length(qcsamples), '/', nrow(pvr), ' samples passed QC.', sep=''))

  passedQC = res[,c(1:5,grep(paste(qcsamples,collapse='|'), colnames(res)))]
  failedQC = res[,unique(c(1:5,grep(paste(qcsamples,collapse='|'), colnames(res), invert=T)) )]
  if (ncol(failedQC) <= 5) failedQC = NULL

  return(list('seg.vals'=passedQC, 'residuals'=pvr, 'prepped.data'=data, 'seg.plots'=plist, 'genome.coverage'=coverage, 'failedQC'=failedQC))
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
#' @param data a data frame of segmented values (from segmentRawData) 
#' @param size tile size to use across the genome, DEF=5e6, options are integer values or "arms"
#' @param build genome build DEF=hg19
#' @return matrix of weighted mean segmented values per tile
#'
#' @author skillcoyne
#' @export
tileSegments<-function(data, size=5e6, build='hg19', verbose=T) {
  require(tibble)

  if (!is.numeric(size) & size != 'arms')
    stop("Size must be numeric, or 'arms'")

  if (!is.tibble(data)) data = as_tibble(data)

  descCols = sort(union(grep('chr|arm|start|end|probes', colnames(data), ignore.case=T), which(!sapply(data, is.numeric))))
  dataCols = c((descCols[length(descCols)]+1):ncol(data))

  chrCol = grep('chr',colnames(data),ignore.case=T, value=T)
  armCol = grep('arm',colnames(data),ignore.case=T,value=T)
  startPos = grep('start',colnames(data),ignore.case=T,value=T)
  endPos = grep('end',colnames(data),ignore.case=T,value=T)

  data = data[which(!data[[chrCol]] %in% c('X','Y')),]
  x1 = data[,c(chrCol, startPos, endPos, colnames(data)[dataCols])]

  tiles = .tile.genome(size, chrInfo(build=build), allosomes=length(which(grepl('X|Y', unique(data[[chrCol]])))) > 0)
  gr = GenomicRanges::makeGRangesFromDataFrame(x1, keep.extra.columns=T, start.field=startPos, end.field=endPos  )

  mergedDf = (do.call(rbind, lapply(tiles, function(tile) {
    cbind('chr'=as.character(seqnames(tile)), as.data.frame(ranges(tile))[1:2])
  }) ))
  mergedDf[colnames(data)[dataCols]] = NA

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
      values = apply(as.data.frame(GenomicRanges::elementMetadata(segments)), 2, weighted.mean, w=weights, na.rm=T)
      mergedDf[rows, names(values)] = values
    }
  }

  message(paste("Mean number of CN segments per genome bin:", round(mean(meanSegs, na.rm=T), 2), "median:", round(median(meanSegs, na.rm=T), 2)))


  rownames(mergedDf) = paste(mergedDf$chr, ':', mergedDf$start, '-', mergedDf$end, sep='')
  mergedDf[,c(1:3)] = NULL
  mergedDf = apply(mergedDf,2, function(x) {
    x[is.na(x)] = median(x,na.rm=T)
    return(x)
  })

  # Not sure if this should be NA or 0
  #mergedDf[is.na(mergedDf)] = 0
  return(t(mergedDf))
}


#' Recommend that this not be used directly. segmentRawData generates these plots as part of the method.
.plotSegmentedGenome<-function(fitted, segmented,window.depths.std,probes.min=NULL) {
  require(ggplot2)

  chr.info = chrInfo(build='hg19')[1:22,]
  chr.info$chrom = sub('chr','', chr.info$chrom)
  chr.info$chrom = factor(chr.info$chrom, levels=c(1:22),ordered = T)

  fitted$chrom = factor(fitted$chrom, levels=levels(chr.info$chrom))
  segmented$chrom = factor(segmented$chrom, levels=levels(chr.info$chrom))

  df = cbind.data.frame('chrom'=fitted$chrom, 'position'=fitted$end, 'seg.cov'=window.depths.std[,1])
  df2 = cbind.data.frame('chrom'=segmented$chrom, 'start'=segmented$start.pos, 'end'=segmented$end.pos, 'seg.val'=segmented[,ncol(segmented)])

  df$chrom = factor(df$chrom, levels=c(1:22), ordered=T)
  df2$chrom = factor(df2$chrom, levels=c(1:22), ordered=T)

  ggplot(chr.info, aes(x=1:chr.length)) +
    ylim(0,4) + facet_grid(~chrom, space='free_x', scales='free_x') +
    geom_point(data=df, aes(x=position, y=seg.cov), color='darkred', alpha=.4) +
    geom_segment(data=df2, aes(x=start, xend=end, y=seg.val, yend=seg.val), color='green3', lwd=5) +
    geom_hline(data=df2, aes(yintercept=median(seg.val)), color='white', linetype='dashed') +
    labs(x='chromosome', y='segmented coverage') + theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0,'lines'))
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
# renamed from prep.pcf.data
.prepRawSWGS<-function(raw.data,fit.data,blacklist, logTransform=F) {
  if (is.null(blacklist) | ncol(blacklist) < 3)
    stop('Blacklisted regions missing or incorrectly formatted.\nExpected columnes: chromosome start end')

  countCols = grep('loc|feat|chr|start|end', colnames(fit.data), invert=T)

  rows = which(fit.data[[1]] %in% raw.data[[1]])
  fit.data = fit.data[rows,]

  window.depths = as.vector(as.matrix(raw.data[,countCols]))/as.vector(as.matrix(fit.data[,countCols]))

  # QDNAseq does this in 'correctBins' but we don't use that method so added here
  negs = which(as.vector(as.matrix(fit.data[,countCols])) <= 0)
  if (length(negs) > 0) window.depths[negs] = 0

  window.depths = matrix(window.depths, ncol=length(countCols))

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

  return(list('data'=data,'sdevs'=sdevs, 'good.bins'=good.bins, 'window.depths.standardised'=window.depths.standardised, 'fit.data'=fit.data))
}
