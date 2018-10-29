
#' Times per sample are determined by the order of the sample as given by the demoFile, or by the order of the samples in the dataset.
#' @name rxRules
#' @return A table describing recommendation rules
#' @author skillcoyne
#' @export
rxRules<-
  tibble(
    'Rule' = c(1:4),
    'Rx' = c('Immediate RFA', 'Recheck in 6-12 months',
             'Recheck in 12-24 months','Regular surveillance in 3-5 years'),
    'Description' = c(
      'HGD or IMC diagnosis or more than one consecutive high risk predictions.',
      'One high risk prediction or an aberrant p53 IHC.',
      'One or more moderate risk predictions.',
      'Two or more consecutive low risk predictions.')
  )


pi.hat<-function(x) exp(x)/(1+exp(x))

#' Main method that should be called to use the trained model.
#' @name predictRisk
#' @param path where QDNAseq raw and corrected read files are REQ
#' @param raw.file.grep for the filename of the raw reads file, default is 'raw.*read'
#' @param corrected.file.grep for the filename of the corrected reads file, default is 'corr|fitted'
#' @return BarrettsRiskRx object. Additional methods that take the object provide information (predictions, rx, segmentedValues, plotSegmentData, adjustRisk, sampleResiduals)
#'
#' @author skillcoyne
#' @export
predictRisk<-function(path='.', raw.file.grep='raw.*read', corrected.file.grep='corr|fitted', cache.dir=NULL, verbose=T) {

  if (is.null(cache.dir))
    cache.dir = getcachedir()

  if (!exists('z.mean')) stop("sysdata file failed to load")
    
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

  if (length(which(sapply(raw.data[,sampleCols], is.integer))) < length(sampleCols))
    stop(paste("Sample columns in raw read count file contain non-integer data."))

  if (length(which(sapply(fit.data[,sampleCols], is.double))) < length(sampleCols))
    warning(paste("Sample columns in fitted read count file contain integer data, fractional values expected."))

  if (verbose) message(paste(length(sampleCols), "sample(s) loaded from",rawFile))

  segmented = segmentRawData(raw.data,fit.data,cache.dir=cache.dir,verbose=verbose)
  psp = predictRiskFromSegments(segmented, verbose=verbose)
  
  return(psp)
}

#' Get the per-sample residuals calculated from the segmentation phase.
#' @name predictRiskFromSegments
#' @param SegmentedSWGS object
#' @return Predict segmented data using included model
#'
#' @author skillcoyne
#' @export
predictRiskFromSegments<-function(swgsObj, verbose=T) {
  if (class(swgsObj)[1] != 'SegmentedSWGS')
    stop("SegmentedSWGS object missing")
  
    mergedDf = tryCatch({
      segtiles = tileSegments(swgsObj, size=5e6,verbose=verbose)
      for (i in 1:ncol(segtiles$tiles))
        segtiles$tiles[,i] = unit.var(segtiles$tiles[,i], z.mean[i], z.sd[i])
      
      armtiles = tileSegments(swgsObj, size='arms',verbose=verbose)
      for (i in 1:ncol(armtiles$tiles))
        armtiles$tiles[,i] = unit.var(armtiles$tiles[,i], z.arms.mean[i], z.arms.sd[i])
      
      cx.score = unit.var(scoreCX(segtiles$tiles,1), mn.cx, sd.cx)
      mergedDf = subtractArms(segtiles$tiles, armtiles$tiles)
      
      cbind(mergedDf, 'cx'=cx.score)
  }, error = function(e) {
    msg = paste("ERROR tiling segmented data:", e)
    print(msg)
  })
  
  sparsed_test_data = Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                             dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
  for(i in colnames(mergedDf)) sparsed_test_data[,i] = mergedDf[,i]
  
  coef.error = .bootstrap.coef.stderr()
  
  Xerr = cbind(segtiles$error, armtiles$error)
  Xerr_diag = apply(Xerr, 1, function(x) { diag(as.matrix(x^2)) })
  # Not sure aobut this step...
  covB = cov(coef.error[,'jack.se',drop=F])[1]
  X = t(mergedDf[,coef.error$coef])
  B = coef.error[,'1',drop=F]
  
  Var_rr = apply(X,2,function(x) sum(x*covB*x))+sapply(Xerr_diag, function(x) sum(B*x*B))
  perSampleError = sqrt( Var_rr )

  RR = predict(fitV, newx=sparsed_test_data, s=lambda, type='link')
  probs = pi.hat(RR)
  
  preds = tibble('Sample'=rownames(probs), 'Probability'=round(probs[,1],2), 
                 'Relative Risk'=RR[,1], 
                 'Risk'=sapply(probs[,1], .risk))
  
  psp = list('predictions'=preds, 'segmented'=swgsObj, 'per.sample.error'=perSampleError, 'tiles'=mergedDf)
  class(psp) <- c('BarrettsRiskRx', class(psp))
  return(psp)
}

#' Calculate upper and lower boudaries for the absolute risk (probabilities) provided by predictRisk
#' @name absoluteRiskCI
#' @param BarrettsRiskRx object
#' @return Data frame of per-sample low and high absolute risk with risk categories applied
#' @author skillcoyne
#' @export
absoluteRiskCI<-function(psp) {
  if (length(which(class(psp) %in% c('BarrettsRiskRx'))) <= 0)
    stop("BarrettsRiskRx object required")
  
  low = round(pi.hat(psp$predictions$`Relative Risk`-psp$per.sample.error),2)
  high = round(pi.hat(psp$predictions$`Relative Risk`+psp$per.sample.error),2)
  preds = tibble('Sample'=psp$predictions$Sample, 
                 'CI.low'=low,
                 'Risk.low'=sapply(low, .risk),
                 'CI.high'=high,
                 'Risk.high'=sapply(high, .risk))
  return(preds)
}

#' Calculate upper and lower boudaries for the relative risk provided by predictRisk
#' @name relativeRiskCI
#' @param BarrettsRiskRx object
#' @return Data frame of per-sample low and high relative risk
#'
#' @author skillcoyne
#' @export
relativeRiskCI<-function(psp) {
  if (length(which(class(psp) %in% c('BarrettsRiskRx'))) <= 0)
    stop("BarrettsRiskRx object required")
  
  low = (psp$predictions$`Relative Risk`-psp$per.sample.error)
  high = (psp$predictions$`Relative Risk`+psp$per.sample.error)
  preds = tibble('Sample'=rownames(probs), 
                 'CI.RR.low'=low,
                 'CI.RR.high'=high)
  return(preds)
}


#' Get the per-sample residuals calculated from the segmentation phase.
#' @name sampleResiduals
#' @param BarrettsRiskRx object
#' @return Data frame of per-sample residual variance, used in QC of samples.
#' @author skillcoyne
#' @export
sampleResiduals<-function(brr) {
  
  if (length(which(class(brr) %in% c('BarrettsRiskRx', 'SegmentedSWGS'))) <= 0)
    stop("BarrettsRiskRx or SegmentedSWGS required")
  if ('SegmentedSWGS' %in% class(brr) )
    return(brr$residuals)
  else 
    return(brr$segmented$residuals)  
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
#' TODO This function presume that each sample is an independent timepoint. As multiple samples are typically collected per timepoint I need to allow the user to indicate that somehow and account for that in the recommendation using the maximum risk.
#' @name rx
#' @param BarrettsRiskRx object REQ
#' @param sample.info Data frame from the 'loadSampleInformation()' function that contains a per-sample entry 
#' @return A table that includes recommendations per timepoint.
#'
#' @author skillcoyne
#' @export
rx<-function(brr, sample.info=NULL) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object required")

  pR = brr$predictions

  if (!is.null(sample.info)) 
    pR = .setUpRxTable(sample.info, pR)
  
  ## TODO check endoscopy number, get max prediction per endo? etc
  

  pR$Risk = factor(pR$Risk, levels=c('Low','Moderate','High'), ordered=T)
  pR$Sample = as.character(pR$Sample)
  p53Col = grep('p53', colnames(pR), value=T, ignore.case=T)
  pathCol = grep('path', colnames(pR), value=T, ignore.case=T)

  rules = as.data.frame(matrix(ncol=3, nrow=nrow(pR)-1, dimnames=list(c(), c('Time 1','Time 2','rule'))))
  # Consecutive
  for (i in 1:(nrow(pR)-1)) {
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

# Current rules
.rule.rx<-function(n) {
  return( subset(rxRules, Rule == n)$Rx )
}


.risk<-function(p) {
  if (!is.numeric(p) | (p > 1 | p < 0) ) stop("Numeric probability between 0-1 required")
  return(as.character(subset(pred.confidence, p <= r2 & p >= r1)$Risk))
}


.readFile<-function(file, ...) {
  if ( tools::file_ext(file) %in% c('xlsx','xls') ) {
    data = readxl::read_xlsx(file,1)
  } else {
    data = readr::read_tsv(file, col_names=T, trim_ws=T, ...)  
  }
  return(data)
}

# Loads a data file with per-sample information on pathology, p53 IHC, etc
.setUpRxTable<-function(sample.info, predDf) {
  if (nrow(sample.info) < nrow(predDf)) {
    warning(paste("Fewer samples in",file,"than in the data files. Ignoring demo data."))
    return(NULL)
  }

  preds = left_join(predDf, sample.info, by='Sample')
    
  return(preds)
}
