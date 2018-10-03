#' Main method that should be called to use the trained model.
#' @name predictRisk
#' @param path where QDNAseq raw and corrected read files are REQ
#' @param raw.file.grep for the filename of the raw reads file, default is 'raw.*read'
#' @param corrected.file.grep for the filename of the corrected reads file, default is 'corr|fitted'
#' @return BarrettsRiskRx object. Additional methods that take the object provide information (predictions, rx, segmentedValues, plotSegmentData, adjustRisk, sampleResiduals)
#'
#' @author skillcoyne
#' @export
predictRisk<-function(path='.', raw.file.grep='raw.*read', corrected.file.grep='corr|fitted', verbose=T) {

  rawFile = grep(raw.file.grep, list.files(path, pattern='txt', full.names=T), value=T, ignore.case=T)
  corrFile = grep(corrected.file.grep, list.files(path, pattern='txt', full.names=T), value=T,ignore.case=T)

  if (length(rawFile) <= 0 | length(corrFile) <= 0)
    stop(paste("Raw reads file or corrected reads file not found in directory:", path))

  if (length(rawFile) > 1 | length(corrFile) > 1)
    stop(paste("Only one raw and fitted file expected, multiple found in: ", path))

  if (verbose)
    message(paste("Raw reads file: ", rawFile, "\n", "Corrected reads file: ", corrFile, "\n", sep=''))

  raw.data = as_tibble(data.table::fread(rawFile,stringsAsFactors=F, showProgress=verbose))
  fit.data = as_tibble(data.table::fread(corrFile,stringsAsFactors=F,showProgress=verbose))

  if (length(which(dim(raw.data) == dim(fit.data))) != 2)
    stop("Raw and corrected read files do not have the same number of rows or columns.")

  descCols = grep('loc|chr|start|end',colnames(raw.data))
  if (length(descCols) != 4)
    stop("File needs the following genome information columns: location chrom start end.\nAll additional column should be samples.")

  sampleCols = grep('loc|chr|start|end',colnames(raw.data), invert=T)
  if (length(sampleCols) < 1)
    stop("No sample columns found after genome information columns.")

  if(length(which(sapply(raw.data[,sampleCols], is.integer))) < length(sampleCols))
    stop(paste("Sample columns in raw read count file contain non-integer data."))

  if(length(which(sapply(fit.data[,sampleCols], is.double))) < length(sampleCols))
    stop(paste("Sample columns in fitted read count file contain integer data, fractional values expected."))

  if (verbose) message(paste(length(sampleCols), "sample(s) loaded from",rawFile))

  segmented = segmentRawData(raw.data,fit.data,verbose=verbose)

  segtiles = tileSegments(segmented$seg.vals, size=5e6,verbose=verbose)
  for (i in 1:ncol(segtiles))
    segtiles[,i] = unit.var(segtiles[,i], z.mean[i], z.sd[i])

  armtiles = tileSegments(segmented$seg.vals, size='arms',verbose=verbose)
  for (i in 1:ncol(armtiles))
    armtiles[,i] = unit.var(armtiles[,i], z.arms.mean[i], z.arms.sd[i])

  cx.score = unit.var(scoreCX(segtiles,1), mn.cx, sd.cx)

  mergedDf = subtractArms(segtiles, armtiles)
  mergedDf = cbind(mergedDf, 'cx'=cx.score)

  sparsed_test_data = Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                             dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
  for(i in colnames(mergedDf)) sparsed_test_data[,i] = mergedDf[,i]

  probs = predict(fitV, newx=sparsed_test_data, s=lambda, type='response')
  RR = predict(fitV, newx=sparsed_test_data, s=lambda, type='link')

  preds = tibble('Sample'=rownames(probs), 'Probability'=probs[,1],
                 'Relative Risk'=RR[,1], 'Risk'=sapply(probs[,1], .risk))

  psp = list('predictions'=preds, 'segmented'=segmented)
  class(psp) <- c('BarrettsRiskRx', class(psp))
  return(psp)
}


#' Get the per-sample residuals calculated from the segmentation phase.
#' @name sampleResiduals
#' @param BarrettsRiskRx object
#' @return Data frame of per-sample residual variance, used in QC of samples.
#'
#' @author skillcoyne
#' @export
sampleResiduals<-function(brr) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object missing")
  brr$segmented$residuals
}


#' Get the samplenames 
#' @name sampleNames
#' @param BarrettsRiskRx object
#' @param passQC logical indicating whether to return the samples that passed (T) or failed (F) QC. DEF=T
#' @return Vector of samplenames
#'
#' @author skillcoyne
#' @export
sampleNames<-function(brr,passQC=T) {
  all = as.character(sampleResiduals(brr)$sample)
  df = brr$segmented$seg.vals  
  if (!passQC) df = brr$segmented$failedQC
  return(intersect(colnames(df), all))
}



#' Get the segmented values for samples from the segmentation phase.
#' @name segmentedValues
#' @param BarrettsRiskRx object
#' @param passQC logical indicating whether to return the samples that passed (T) or failed (F) QC. DEF=T
#' @return Data frame of segmented values for samples that passed QC
#'
#' @author skillcoyne
#' @export
segmentedValues<-function(brr, passQC=T) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object missing")
  df = brr$segmented$seg.vals  
  if (!passQC) df = brr$segmented$failedQC
  return(df)
}

#' Plot genome-wide raw and segmented values for all samples in the segmentation phase.
#' @name plotSegmentData
#' @param BarrettsRiskRx object
#' @return ggplot of raw and segmented values
#'
#' @author skillcoyne
#' @export
plotSegmentData<-function(brr) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object missing")
  do.call(gridExtra::grid.arrange, c(brr$segmented$seg.plots, ncol=1))
}

#' Use with caution. These are based on estimates of what we think the risk might really be #' and adjusting the resulting risk based on those estimates.
#' @name adjustRisk
#' @param BarrettsRiskRx object
#' @param offset c(min,mean,max)
#' @return BarrettsRiskRx object with adjusted predictions and recommentations.
#'
#' @author skillcoyne
#' @export
adjustRisk <- function(brr, offset=c('min','mean','max')) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object missing")

  if (length(grep('segmented',names(brr))) <= 0)
    stop("Original BarrettsRiskRx object required. Adjusted object should not be re-adjusted.")

  total = cbind.data.frame('NP'=43,'P'=45)
  cases = total[['P']]
  
  mn = round((cases/(0.01*100))*100)
  m = round((cases/(0.0225*100))*100)
  mx = round((cases/(0.035*100))*100)
  
  offset = switch(offset,
         min=log(cases/mn),
         max=log(cases/mx),
         mean=log(cases/m))

  RR = brr$predictions$`Relative Risk`

  adjustedRiskRx = brr$predictions

  adjustedRiskRx$Probability = 1/(1+exp(-RR+abs(offset)))
  adjustedRiskRx$`Relative Risk` = RR+offset
  adjustedRiskRx$Risk = sapply(adjustedRiskRx$Probability, .risk)

  psp = list('predictions'=adjustedRiskRx)
  class(psp) <- c('BarrettsRiskRx', class(psp))
  return(psp)
}

#' Get predictions from the model
#' @name predictions
#' @param BarrettsRiskRx object REQ
#' @return predictions data frame
#'
#' @author skillcoyne
#' @export
predictions<-function(brr) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object missing")
  brr$predictions
}

#' Times per sample are determined by the order of the sample as given by the demoFile, or by the order of the samples in the dataset.
#' @name rx
#' @param BarrettsRiskRx object REQ
#' @param demoFile for the tab-delimited file that contains a per-sample entry for p53 IHC, Barrett's segment length, patient gender OPTIONAL
#' @return A table that includes recommendations per timepoint.
#'
#' @author skillcoyne
#' @export
rx<-function(brr, demoFile=NULL) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object required")

  pR = brr$predictions

  if (!is.null(demoFile))
    pR = .loadDemoData(demoFile, pR$Sample, pR)

  pR$Risk = factor(pR$Risk, levels=c('Low','Moderate','High'), ordered=T)
  pR$Sample = as.character(pR$Sample)
  p53Col = grep('p53', colnames(pR), value=T, ignore.case=T)
  pathCol = grep('path', colnames(pR), value=T, ignore.case=T)

  rules = as.data.frame(matrix(ncol=3, nrow=nrow(pR), dimnames=list(c(), c('Time 1','Time 2','rule'))))
  # Consecutive
  for (i in 1:nrow(pR)) {
    risks = table(pR$Risk[i:(i+1)])
    p53 = NULL
    if (length(p53Col) > 0) {
      pR[[p53Col]] = factor(pR[[p53Col]], levels=c(0,1))
      p53 = table(pR[[p53Col]][i:(i+1)])
    }

    rule = 'None'
    if ( risks['High'] == 2 || (length(pathCol) > 0 && length(which(grepl('HGD|IMC', pR[[pathCol]][i:(i+1)]))) > 0) ) {
      rule = 1
    } else if ( risks['High'] == 1 || (!is.null(p53) && p53['1'] > 0) ) {
      rule = 2
    } else if ( risks['Moderate'] > 0 || (risks['Low'] ==1 && nrow(pR) == 1)) {
      rule = 3
    } else if ( risks['Low'] == 2 ) {
      rule = 4
    }
    rules[i,] = cbind( pR$Sample[i], pR$Sample[(i+1)], as.integer(rule) )
  }
  rules$rule = as.integer(rules$rule)
  rules$Rx = sapply(rules$rule, .rule.rx)

  return(rules)
}


#' Model predictions plot
#' @name plotModelPredictions
#' @param type RR (relative risk) or P (probability, DEF)
#' @return ggplot object 
#'
#' @author skillcoyne
#' @export
plotModelPredictions<-function(type='P') {
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

# Current rules
.rule.rx<-function(n) {
  rr = c('Immediate RFA', 'Recheck 6-12 months',
         'Recheck 12-24 months','Regular surveillance 3-5 years')
  return(rr[n])
}


.risk<-function(p) {
  if (!is.numeric(p) | (p > 1 | p < 0) ) stop("Numeric probability between 0-1 required")
  return(as.character(subset(pred.confidence, p <= r2 & p >= r1)$Risk))
}

# Loads a data file with per-sample information on pathology, p53 IHC, etc
.loadDemoData<-function(file, samplenames, predDf) {
  demo = as.data.frame(data.table::fread(file))

  if (nrow(demo) > length(samplenames)) {
    warning(paste("Not all samples in",file,"are represented in the data files."))
    demo = subset(demo, Sample %in% samplenames)
  }

  if (nrow(demo) < length(samplenames)) {
    warning(paste("Fewer samples in",file,"than in the data files. Ignoring demo data."))
    return(NULL)
  }

  path = grep('path',colnames(demo),ignore.case=T,value=T)
  p53 = grep('p53',colnames(demo),ignore.case=T,value=T)

  if( length(path) <= 0 | length(p53) <= 0 ) {
    warning(paste("p53 IHC or pathology not available, recommendations will be limited to available information and risk prediction."))
  }

  beLen = grep('length',colnames(demo),ignore.case=T,value=T)
  gender = grep('gender|sex',colnames(demo),ignore.case=T,value=T)

  preds = base::merge(predDf, demo[,c('Sample',path,p53,beLen,gender)],by='Sample')

  #intersect(demo$Sample, preds$Sample)

  return(preds)
}
