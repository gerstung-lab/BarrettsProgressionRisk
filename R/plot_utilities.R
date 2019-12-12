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
#' @param as Return type, plot or list of plots
#' @return ggplot of raw and segmented values
#'
#' @author skillcoyne
#' @export
plotSegmentData<-function(brr, as=c('plot','list')) {
  if (length(which(class(brr) %in% c('BarrettsRiskRx', 'SegmentedSWGS'))) <= 0)
    stop("BarrettsRiskRx or SegmentedSWGS required")
  rettype = match.arg(as)  
  
  if ('SegmentedSWGS' %in% class(brr))  {
    plotlist = brr$seg.plots
  } else {
    plotlist = brr$segmented$seg.plots  
  }

  if (rettype == 'plot')
    return(do.call(gridExtra::grid.arrange, c(plotlist, ncol=1)))
  else 
    return(plotlist)
}


#' Plot genome-wide coverage from adjusted raw data
#' @name plotCorrectedCoverage
#' @param BarrettsRiskRx or SegmentedSWGS object
#' @return ggplot of raw and segmented values
#'
#' @author skillcoyne
#' @export
plotCorrectedCoverage<-function(brr, as=c('plot','list')) {
  if (length(which(class(brr) %in% c('BarrettsRiskRx', 'SegmentedSWGS'))) <= 0)
    stop("BarrettsRiskRx or SegmentedSWGS required")
  rettype = match.arg(as)  
  
  if ('SegmentedSWGS' %in% class(brr))  {
    plotlist = brr$cv.plot
  } else {
    plotlist = brr$segmented$cv.plot  
  }
  
  if (rettype == 'plot')
    return(do.call(gridExtra::grid.arrange, c(plotlist, ncol=1)))
  else 
    return(plotlist)
}

#' Predictions risk calibration plot
#' @name showPredictionCalibration
#' @return ggplot object 
#'
#' @author skillcoyne
#' @export
showPredictionCalibration<-function(df=NULL) {
  if (is.null(df)) df = BarrettsProgressionRisk:::be_model$pred.confidence
  
  mm = range(df[c('r1','r2')])
  
  cuts = seq(mm[1], mm[2], by=df$r1[2])

  plot.theme = theme(text=element_text(size=12), panel.background=element_blank(), strip.background =element_rect(fill="white"),  
                     strip.text = element_text(size=12), 
                     axis.line=element_line(color='black'), panel.grid.major=element_line(color='grey90'),
                     panel.border = element_rect(color="grey", fill=NA, size=0.5), panel.spacing = unit(0.1, 'lines')  ) 
  
  ggplot(df, aes(mn, perc)) + 
    geom_rect(aes(xmin=r1, xmax=r2, ymin=0,ymax=1, fill=Risk), alpha=0.6) + 
    scale_fill_manual(values=riskColors(), limits=levels(df$Risk) ) +
    geom_vline(xintercept=cuts[2:(length(cuts)-1)], color='grey88') +
    geom_smooth(method='lm',formula=y~x, color='grey39', linetype='dashed', fill='grey88', size=0.5, fullrange=T) + 
    geom_point() + geom_errorbar(aes(ymin=ci.low, ymax=ci.high), size=0.5, width=0.01) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    scale_x_continuous(expand=c(0,0),limits=c(-0.5,1.5), breaks=cuts, labels=cuts) + 
    scale_y_continuous(expand=c(0,0),limits=c(-0.5,1.5), breaks=cuts, labels=cuts) +
    plot.theme + theme(legend.position = 'bottom') + labs(x='mean(Absolute Risk)', y='Progressor:Non-Progressor', title='Risk Calibration') 
}

#' Get the colors used for risk classes 'low', 'moderate','high'
#' @name riskColors
#' @return tibble
#'
#' @author skillcoyne
#' @export
riskColors<-function() {
  rc = tibble(Low = "#4575B4", Moderate="#FEE090",High="#A50026")
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

#' Per sample over time tile risk plot
#' @name patientRiskTilesPlot
#' @param brr BarrettsRiskRx object or tibble with Risk column
#' @param col time column name, default 'Endoscopy'
#' @param direction to plot the time, fwd or rev  (forward or reverse)
#' @return ggplot object
#'
#' @author skillcoyne
#' @export
patientRiskTilesPlot<-function(brr, col='Endoscopy', direction=c('fwd','rev')) {
  dir = match.arg(direction)  
  
  if (length(which(class(brr) %in% c('BarrettsRiskRx'))) > 0) {
    preds = brr$per.sample
  } else {
    preds = brr
    if (length(which(colnames(preds) %in% c(col, 'Risk'))) < 2) stop("BarrettsRiskRx object required, or tibble with numeric/date and Risk (Low,Moderate,High) columns")
  }

  if ('GEJ.Distance' %in% colnames(preds) & !is.factor(preds$GEJ.Distance)) {
    preds$GEJ.Distance = fct_rev(factor(preds$GEJ.Distance, ordered=T))
  } else if (!'GEJ.Distance' %in% colnames(preds)) {
    preds$GEJ.Distance = 1
  }
  preds = preds %>% mutate_if(is.numeric, list(~factor(.,ordered=T)))
  
  p = ggplot(preds, aes_string(col, 'GEJ.Distance')) +
    geom_tile(aes(fill=Risk), color='white',size=2) + 
    scale_fill_manual(values=riskColors(), limits=names(riskColors())) +
    labs(y='Esophageal Location (GEJ...)')
  
  if (dir == 'rev') p = p + scale_x_discrete(limits=rev(levels(preds[[col]])))

  if ('Pathology' %in% colnames(preds)) {
    p = p + geom_point(aes(shape=Pathology), fill='white', color='white', size=8) + 
        scale_shape_manual(values=c(1,0,15,24,25), limits=c('NDBE','ID','LGD','HGD','IMC'), labels=c('NDBE','ID','LGD','HGD','IMC'), guide=guide_legend(override.aes=list(fill='white', color='white')))
    }
  
  if ('Patient' %in% colnames(preds)) p = p + labs(title=unique(preds$Patient))
  
  p + theme_bw() + theme(legend.key=element_rect(fill='grey39'), panel.background=element_rect(colour = 'black'), panel.grid.major=element_blank(), panel.spacing = unit(0.2, 'lines'), panel.border = element_rect(color="black", fill=NA, size=0.5), legend.position = 'bottom'  ) 
}


#' Over time per endoscopy risk plot
#' @name patientEndoscopyPlot
#' @param type BarrettsRiskRx object
#' @return ggplot object
#'
#' @author skillcoyne
#' @export
patientEndoscopyPlot<-function(brr) {
  if (length(which(class(brr) %in% c('BarrettsRiskRx'))) <= 0)
    stop("BarrettsRiskRx required")
  
  preds = absoluteRiskCI(brr)
  preds = preds %>% rowwise() %>% dplyr::mutate( img=printRisk(Probability*100,CI.low*100,CI.high*100,Risk) )
  
  ggplot(preds, aes(Endoscopy, Probability)) + ylim(0,1) +
    geom_line(color='grey') + 
    geom_errorbar(aes(ymin=CI.low,ymax=CI.high, color=Risk), width=5, show.legend=F) + 
    geom_point(aes(color=Risk), size=5) + 
    scale_color_manual(values=riskColors(), limits=names(riskColors())) + 
    labs(y='Absolute Risk', x='Endoscopy Date',title='Absolute risks over time') + theme_bw() + theme(legend.position='bottom')
}


#' Windowed, scaled CN values across the genome
#' @name copyNumberMountainPlot
#' @param type BarrettsRiskRx object
#' @param annotate If true, only segments in the non-zero coefficients are highlighted (DEF=T)
#' @param as plot or list of plots
#' @return list of ggplot objects per sample
#'
#' @author skillcoyne
#' @export
copyNumberMountainPlot<-function(brr,annotate=T, legend=T,  as=c('plot','list')) {
  if (length(which(class(brr) %in% c('BarrettsRiskRx'))) <= 0)
    stop("BarrettsRiskRx required")
  rettype = match.arg(as)
  
  mp = BarrettsProgressionRisk:::mountainPlots(brr$tiles, as.matrix(coef(brr$be.model$fit, brr$be.model$lambda)), brr$be.model$cvRR, brr$segmented$chr.build.info, annotate=T) 
  ht = length(mp$plot.list)*2
  
  if (rettype == 'list') {
    plist = mp$plot.list
    if (legend)
      plist = lapply(mp$plot.list, function(p) gridExtra::arrangeGrob(p, mp$legend,  heights=c(ht = length(mp$plot.list)*2,1), ncol=1))
    return( plist )
  } else {
    p = do.call(gridExtra::arrangeGrob, (mp$plot.list))
    
    if (legend)
      return(gridExtra::grid.arrange(p,gridExtra::arrangeGrob(mp$legend), heights=c(ht = length(mp$plot.list)*2,1), ncol=1))
    else
      return(gridExtra::grid.arrange(p, ncol=1))
  }
}

mountainPlots<-function(tiles, coefs, cvRR, build, annotate=T) {
  pal = c('#238B45', 'grey','#6A51A3')
  chr.info = BarrettsProgressionRisk:::chrInfo(build=build)
  
  samples = rownames(tiles)
  locs = BarrettsProgressionRisk:::get.loc(tiles[,-ncol(tiles),drop=F])

  #coefs =  as.matrix(coef(brr$be.fit$fit, brr$be.fit$lambda))
  if (!is.null(coefs)) {
    coefs = coefs[which(coefs != 0),][-1]
    
    coefs = tibble::enframe(coefs, name='feature') %>% dplyr::rename(coef = 'value') %>%
      separate(feature, c('chr','start','end'), sep = ':|-', remove = F) %>% 
      mutate_at(vars(start, end), as.double) %>%       
      mutate(chr = factor(chr, levels=levels(locs$chr), ordered=T))
    
    #if (is.data.frame(cvRR) | is.matrix(cvRR)) cvRR = cvRR[,'cvRR']
  
    cvdf = cvRR %>% separate(label, c('chr','start','end'), sep=':|-', remove=F) %>%
      mutate_at(vars(start, end), as.double) %>%       
      mutate(chr = factor(chr, levels=levels(locs$chr), ordered=T)) %>%
      arrange(chr,start)

    cvdf = full_join(cvdf, coefs, by=c('label'='feature','chr','start','end', 'coef'))
  }
  
  plist = list()
  for (sample in samples) {
    df = bind_cols(locs,tibble('sample'=tiles[sample, -ncol(tiles)]))

    melted = as_tibble(reshape2::melt(id.vars=c('chr','start','end'), df))
    melted = left_join(melted,chr.info[,c('chr','chr.length')],by='chr')
    if (!is.null(coefs)) {
      melted = as_tibble(left_join(melted,cvdf, by=c('chr','start','end')) )
    }

    arms = melted %>% filter(end-start > median(melted$end-melted$start)*2)
    segs = melted %>% filter(end-start <= median(melted$end-melted$start)*2)
    
    cn<-function(value,range) {
      x = 'norm'
      if (value < range[1]) {
        x = 'loss' 
      } else if (value > range[2]) {
        x = 'gain'
      }
      return(x)
    }
  
    segs = segs %>% rowwise() %>% dplyr::mutate( 
      CN = cn(value,c(-1,1)),
      annotate = ifelse(!is.na(coef),cn(value,c(0,0)),'norm'),
    ) %>% mutate(
      CN = factor(CN, levels = c('loss','norm','gain')),
      annotate = factor(annotate, levels = c('loss','norm','gain'))
    )

    arms = arms %>% rowwise() %>% dplyr::mutate( 
        CN = cn(value,c(-1,1)), 
        annotate = ifelse(!is.na(coef),cn(value,c(0,0)),'norm')) %>% 
      filter(!is.na(coef)) %>% mutate(
        CN = factor(CN, levels = c('loss','norm','gain')),
        annotate = factor(annotate, levels = c('loss','norm','gain'))
      )

    p = ggplot(segs, aes(x=chr.length)) + facet_grid(~chr, scales='free_x', space='free_x') 

    if (!annotate) {
      p = p + geom_rect(aes(xmin=1,xmax=chr.length,ymin=-1,ymax=1),fill='grey88',alpha=0.03) + 
        geom_rect(aes(xmin=start,xmax=end,ymin=0, ymax=value, fill=CN)) + 
        scale_fill_manual(values=pal, limits=c('loss','norm','gain'), labels=c('Loss','Normal','Gain'), name='Relative CNA')
    } else {
      p = p + geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=value, fill=annotate)) + 
        geom_rect(data=arms,aes(xmin=start,xmax=end,ymin=0,ymax=value,fill=annotate),alpha=0.5) +
        geom_rect(aes(xmin=1,xmax=chr.length,ymin=-1,ymax=1),fill='grey88',alpha=0.03) +
        scale_fill_manual(values=pal, limits=c('loss','norm','gain'), labels=c('Loss','Normal','Gain'), name='Model Features')
    }
    p = p + labs(x='Chromosomes', y='Relative CNA',title=sample) +
      theme_bw() + theme(axis.text.x=element_blank(), panel.spacing.x=unit(0,'lines'), panel.grid=element_blank(), panel.border=element_rect(linetype='solid', color='grey39', fill=NA), legend.position = 'bottom' ) 
    
    legend = .get.legend(p)
    
    plist[[sample]] = p + theme(legend.position = 'none')
  }
  legend = gridExtra::arrangeGrob(legend)
  
  return(list('plot.list'=plist, 'legend'=legend))
}


.get.legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



