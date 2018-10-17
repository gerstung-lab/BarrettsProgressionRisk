require(ggplot2)
require(RColorBrewer)




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


#' Plot genome-wide raw and segmented values for all samples in the segmentation phase.
#' @name plotSegmentData
#' @param BarrettsRiskRx or SegmentedSWGS object
#' @return ggplot of raw and segmented values
#'
#' @author skillcoyne
#' @export
plotSegmentData<-function(brr) {
  if (length(which(class(brr) %in% c('BarrettsRiskRx', 'SegmentedSWGS'))) <= 0)
    stop("BarrettsRiskRx or SegmentedSWGS required")
  if ('SegmentedSWGS' %in% class(brr) )
    do.call(gridExtra::grid.arrange, c(brr$seg.plots, ncol=1))
  else 
    do.call(gridExtra::grid.arrange, c(brr$segmented$seg.plots, ncol=1))
}


#' Plot genome-wide coverage from adjusted raw data
#' @name plotCorrectedCoverage
#' @param BarrettsRiskRx or SegmentedSWGS object
#' @return ggplot of raw and segmented values
#'
#' @author skillcoyne
#' @export
plotCorrectedCoverage<-function(brr) {
  if (length(which(class(brr) %in% c('BarrettsRiskRx', 'SegmentedSWGS'))) <= 0)
    stop("BarrettsRiskRx or SegmentedSWGS required")
  if ('SegmentedSWGS' %in% class(brr) )
    return(brr$cv.plot)
  else 
    return(brr$segmented$cv.plot)
}

#' Predictions risk calibration plot
#' @name showPredictionCalibration
#' @return ggplot object 
#'
#' @author skillcoyne
#' @export
showPredictionCalibration<-function() {
  plot.theme = theme(text=element_text(size=12), panel.background=element_blank(), strip.background =element_rect(fill="white"),  
                     strip.text = element_text(size=12), 
                     axis.line=element_line(color='black'), panel.grid.major=element_line(color='grey90'),
                     panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing = unit(0.1, 'lines')  ) 
  
  ggplot(pred.confidence, aes(mn, perc)) + geom_smooth(method='lm',formula=y~x, color='grey39', linetype='dashed', size=0.5, fullrange=T) + 
    geom_rect(aes(xmin=r1, xmax=r2, ymin=0,ymax=1, fill=Risk), alpha=0.6) + 
    scale_fill_manual(values=risk.colors, limits=levels(pred.confidence$Risk) ) +
    geom_point() + 
    geom_errorbar(aes(ymin=ci.low, ymax=ci.high), size=0.5, width=0.01) +
    scale_color_manual(values=risk.colors,limits=levels(pred.confidence$Risk) ) + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    scale_x_continuous(expand=c(0,0),limits=c(-0.5,1.5), breaks=cuts, labels=cuts) + 
    scale_y_continuous(expand=c(0,0),limits=c(-0.5,1.5), breaks=cuts, labels=cuts) +
    plot.theme + theme(legend.position = 'bottom') + labs(x='P(Progression)', y='Progressor:Non-Progressor', title='Risk Calibration') 
}

#' Get the colors used for risks
#' @name riskColors
#' @return vector
#'
#' @author skillcoyne
#' @export
riskColors<-function() {
  rc = risk.colors
  names(rc) = c('Low','Moderate','High')
  return(rc)
}


#' Model predictions plot
#' @name showModelPredictions
#' @param type RR (relative risk) or P (probability, DEF)
#' @return ggplot object 
#'
#' @author skillcoyne
#' @export
showModelPredictions<-function(type='P') {
  plot.theme = theme(text=element_text(size=12), panel.background=element_blank(), strip.background =element_rect(fill="white"),  
                     strip.text = element_text(size=12), 
                     axis.line=element_line(color='black'), panel.grid.major=element_line(color='grey90'),
                     panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing = unit(0.1, 'lines')  ) 
  
  p = ggplot(cxPredictions, aes(Prediction)) + geom_histogram(aes(fill=..x..), breaks=cuts, show.legend = F) +
    scale_fill_gradientn(colors = hist.pal,  name='') + 
    plot.theme + labs(title='All samples model predictions', y='n Samples', x='Probability') 
  if (type == 'RR') {
    p = ggplot(cxPredictions, aes(Relative.Risk)) + geom_histogram(aes(fill=..x..), bins=20, show.legend = F) +
      scale_fill_gradientn(colors = hist.pal,  name='') + 
      plot.theme + labs(title='All samples model predictions', y='n Samples', x='Relative Risk') 
  }
  return(p)
}
