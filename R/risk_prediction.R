pi.hat<-function(x) exp(x)/(1+exp(x))


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


#' Main method that should be called to use the trained model.
#' @name predictRisk
#' @info SampleInformation object
#' @param path where QDNAseq raw and corrected read files are REQ
#' @param raw.file.grep for the filename of the raw reads file, default is 'raw.*read'
#' @param corrected.file.grep for the filename of the corrected reads file, default is 'corr|fitted'
#' @return BarrettsRiskRx object. Additional methods that take the object provide information (predictions, rx, segmentedValues, plotSegmentData, adjustRisk, sampleResiduals)
#'
#' @author skillcoyne
#' @export
predictRisk<-function(info, path, raw.file.grep='raw.*read', corrected.file.grep='corr|fitted', cache.dir=NULL, verbose=T) {
  if (!'SampleInformation' %in% class(info)) 
    stop("SampleInformation object from loadSampleInformation(...) required.")

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

  raw.data = readr::read_tsv(rawFile, col_names=T, col_types = cols('chrom'=col_character()))
  fit.data = readr::read_tsv(corrFile, col_names=T, col_type = cols('chrom'=col_character()))

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

  segmented = segmentRawData(info, raw.data,fit.data,cache.dir=cache.dir,verbose=verbose)
  psp = predictRiskFromSegments(segmented, verbose=verbose)
  
  return(psp)
}

#' Get the per-sample residuals calculated from the segmentation phase.
#' @name predictRiskFromSegments
#' @param SegmentedSWGS object
#' @param glmnet model for prediction (Default uses internal model)
#' @param s lambda value for prediction (Default uses internal model)
#' @param size Segment size (Default 5e6)
#' @param tile.mean
#' @param tile.sd
#' @param arms.mean
#' @param arms.sd
#' @param cx.mean
#' @param sd.cx
#' @return Predict segmented data using included model
#'
#' @author skillcoyne
#' @export
predictRiskFromSegments<-function(obj, model=fitV, s=lambda, tile.size=5e6, tile.mean=z.mean, arms.mean=z.arms.mean, tile.sd=z.sd, arms.sd=z.arms.sd, cx.mean=mn.cx, cx.sd=sd.cx, verbose=T) {
  if (class(obj)[1] != 'SegmentedSWGS')
    stop("SegmentedSWGS object missing")
  
    if (nrow(sampleResiduals(obj)) == nrow(sampleResiduals(obj) %>% dplyr::filter(!Pass))) {
      warning('No samples passed QC, no predictions can be made.')
      return(NULL)
    }
  
  if (verbose) 
    message( paste('Using internal glmnet model: ', fitV$nobs == model$nobs)  )
  
  if (fitV$nobs != model$nobs & length(which(z.mean == tile.mean) != length(z.mean))) 
    stop('Using external glmnet model. Tile and arm mean/sd values required to correctly mean adjust tiled data.')


    # Tile, scale, then merge segmented values into 5Mb and arm-length windows across the genome.
    mergedDf = tryCatch({
      segtiles = tileSegments(obj, size=tile.size,verbose=verbose)
      for (i in 1:ncol(segtiles$tiles))
        segtiles$tiles[,i] = unit.var(segtiles$tiles[,i], tile.mean[i], tile.sd[i])
      
      armtiles = tileSegments(obj, size='arms',verbose=verbose)
      for (i in 1:ncol(armtiles$tiles))
        armtiles$tiles[,i] = unit.var(armtiles$tiles[,i], arms.mean[i], arms.sd[i])
      
      cx.score = unit.var(scoreCX(segtiles$tiles,1), cx.mean, cx.sd)
      mergedDf = subtractArms(segtiles$tiles, armtiles$tiles)
      
      cbind(mergedDf, 'cx'=cx.score)
  }, error = function(e) {
    msg = paste("ERROR tiling segmented data:", e)
    stop(msg)
  })
  
  sparsed_test_data = Matrix(data=0, nrow=nrow(mergedDf),  ncol=ncol(mergedDf),
                             dimnames=list(rownames(mergedDf),colnames(mergedDf)), sparse=T)
  for(i in colnames(mergedDf)) sparsed_test_data[,i] = mergedDf[,i]
  
  # get the bootstrap errors for coefficients
  coef.error = .bootstrap.coef.stderr()
  
  # Errors are the per window, weighted mean error of all segments in the bin
  Xerr = cbind(segtiles$error, armtiles$error)
  #if (nrow(Xerr) > 1) {
  Xerr_diag = apply(Xerr, 1, function(x) { diag(as.matrix(x^2)) })
  #} else {
  #  Xerr_diag = diag(as.matrix(Xerr)^2)  
  #}

  # Not sure aobut this step...
  covB = cov(coef.error[,'jack.se',drop=F])[1]
  X = t(mergedDf[,coef.error$coef,drop=F])
  B = coef.error[,'1',drop=F]
  
  Var_rr = apply(X,2,function(x) sum(x*covB*x)) + sapply(Xerr_diag, function(x) sum(B*x*B))
  perSampleError = sqrt( Var_rr )
  perSampleError = tibble('Sample'=names(perSampleError),'Error'=perSampleError)
  
  perSampleError = left_join(perSampleError, obj$sample.info, by='Sample') %>% dplyr::select('Sample','Error','Endoscopy')
  
  # Predict and generate absolute probabilities
  RR = predict(model, newx=sparsed_test_data, s=s, type='link')
  probs = pi.hat(RR)
  
  per.sample.preds = full_join(tibble('Sample'=rownames(probs), 
                                      'Probability'=round(probs[,1],2), 
                                      'Relative Risk'=RR[,1],
                                      'Risk'=sapply(probs[,1], .risk)), 
                               obj$sample.info %>% dplyr::filter(Sample %in% as.character(sampleResiduals(obj) %>% dplyr::filter(Pass) %>% dplyr::select(sample) %>% pull) ), 
                               by='Sample')
  
  per.endo.preds = .setUpRxTablePerEndo(per.sample.preds, 'max', verbose) %>% rowwise() %>% mutate( Risk=.risk(Probability))
  perEndoError = .setUpRxTablePerEndo(perSampleError, 'max', F) %>% dplyr::select('Endoscopy','Error')
  
  psp = list('per.endo'=per.endo.preds, 'per.sample'=per.sample.preds, 'segmented'=obj, 'per.sample.error'=perSampleError, 'per.endo.error'=perEndoError, 'tiles'=mergedDf)
  class(psp) <- c('BarrettsRiskRx', class(psp))
  return(psp)
}

#' Calculate upper and lower boudaries for the absolute risk (probabilities) provided by predictRisk
#' @name absoluteRiskCI
#' @param BarrettsRiskRx object
#' @return Data frame of per-sample low and high absolute risk with risk categories applied
#' @author skillcoyne
#' @export
absoluteRiskCI<-function(psp, by=c('endoscopy','sample'), verbose=T) {
  if (length(which(class(psp) %in% c('BarrettsRiskRx'))) <= 0)
    stop("BarrettsRiskRx object required")
  riskBy = .titleCase(match.arg(by))
  
  if (verbose)
    message(paste('Predictions, risks, and CI per',tolower(riskBy)))

  preds = switch(riskBy,
                 'Sample'=full_join(psp$per.sample, psp$per.sample.error, by=c('Sample','Endoscopy')),
                 'Endoscopy'=full_join(psp$per.endo, psp$per.endo.error, by='Endoscopy'))
  
  low = round(pi.hat(preds$`Relative Risk`-preds$Error),2)
  high = round(pi.hat(preds$`Relative Risk`+preds$Error),2)
  
  preds = add_column(preds, 
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
relativeRiskCI<-function(psp, by=c('endoscopy','sample'), verbose=T) {
  if (length(which(class(psp) %in% c('BarrettsRiskRx'))) <= 0)
    stop("BarrettsRiskRx object required")
  riskBy = .titleCase(match.arg(by))
  
  if (verbose)
    message(paste('Predictions, relative risks, and CI per',tolower(riskBy)))
  
  preds = full_join(psp$per.endo, psp$per.endo.error, by='Endoscopy')
  if (riskBy == 'Sample') 
    preds = full_join(psp$per.sample, psp$per.sample.error, by=c('Sample', 'Endoscopy'))

  low = (preds$`Relative Risk`-preds$Error)
  high = (preds$`Relative Risk`+preds$Error)
  
  preds = add_column(preds, 'CI.RR.low'=low,'CI.RR.high'=high)

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

#' Use with caution! These are based on estimates of what we think the risk might really be and adjusting the resulting risk based on those estimates.
#' @name adjustRisk
#' @param BarrettsRiskRx object
#' @param offset c(min,mean,max)
#' @return list with adjusted predictions and recommentations.
#'
#' @author skillcoyne
#' @export
adjustRisk <- function(brr, offset=c('mean','max','min'), by=c('endoscopy','sample')) {
  offset = match.arg(offset)
  riskBy = .titleCase(match.arg(by)) 
    
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

  preds = switch(riskBy,
                 'Endoscopy'=brr$per.endo,
                 'Sample'=brr$per.sample)
                   
  RR = preds$`Relative Risk`

  adjustedRiskRx = preds

  adjustedRiskRx$Probability = 1/(1+exp(-RR+abs(offset)))
  adjustedRiskRx$`Relative Risk` = RR+offset
  adjustedRiskRx$Risk = sapply(adjustedRiskRx$Probability, .risk)

  recommendations = .apply.rules(adjustedRiskRx,riskBy)
  
  return(list('adj.predictions'=adjustedRiskRx, 'adj.recommendations'=recommendations))
}

#' Get predictions from the model
#' @name predictions
#' @param BarrettsRiskRx object REQ
#' @param c(sample, endo)
#' @return predictions data frame
#'
#' @author skillcoyne
#' @export
predictions<-function(brr, func=c('sample','endoscopy')) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object missing")
  func = match.arg(func)
  
  if (grepl(func, 'sample'))
    return(brr$per.sample)
  if (grepl(func, 'endoscopy'))
    return(brr$per.endo)
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
rx<-function(brr, by=c('endoscopy','sample')) {
  if (class(brr)[1] != 'BarrettsRiskRx')
    stop("BarrettsRiskRx object required")

  riskBy = .titleCase(match.arg(by))

  preds = switch(riskBy,
                 'Endoscopy'=brr$per.endo,
                 'Sample'=brr$per.sample)

  rules = .apply.rules(preds,riskBy)

  return(rules)
}

.apply.rules<-function(preds,riskBy) {
  preds = preds %>% rowwise() %>% dplyr::mutate(Risk = .risk(Probability))
  preds$Risk = factor(preds$Risk, levels=c('Low','Moderate','High'), ordered=T)
  
  p53Col = grep('p53', colnames(preds), value=T, ignore.case=T)
  pathCol = grep('path', colnames(preds), value=T, ignore.case=T)
  
  # Consecutive after sorting by the selected column
  preds = preds %>% arrange(preds[[riskBy]])
  
  rules = as_tibble(matrix(ncol=3, nrow=0, dimnames=list(c(), c('Time 1','Time 2','Rule'))))
  for (i in 1:nrow(preds)) {
    risks = table(preds$Risk[i:(i+1)])
    p53 = NULL
    if (length(p53Col) > 0) {
      preds[[p53Col]] = factor(preds[[p53Col]], levels=c(0,1))
      p53 = table(preds[[p53Col]][i:(i+1)])
    }
    
    rule = 'None'
    if ( risks['High'] == 2 || (length(pathCol) > 0 && length(which(grepl('HGD|IMC', preds[[pathCol]][i:(i+1)]))) > 0) ) {
      rule = 1
    } else if ( risks['High'] == 1 || (!is.null(p53) && p53['1'] > 0) ) {
      rule = 2
    } else if ( risks['Moderate'] > 0 || (risks['Low'] ==1 && nrow(preds) == 1)) {
      rule = 3
    } else if ( risks['Low'] == 2 ) {
      rule = 4
    }
    
    rules = add_row(rules, 'Time 1'=preds[[riskBy]][i], 'Time 2'=preds[[riskBy]][(i+1)], 'Rule'=as.integer(rule)  )
    if (i == nrow(preds)) break;
  }
  rules$Rx = sapply(rules$Rule, .rule.rx)
  
  return(rules)
}

# Current rules
.rule.rx<-function(n) {
  return( subset(rxRules, Rule == n)$Rx )
}


.risk<-function(p) {
  if (!is.numeric(p) | (p > 1 | p < 0) ) stop("Numeric probability between 0-1 required")
  return(as.character(max(subset(pred.confidence, p <= r2 & p >= r1)$Risk)))
}


.readFile<-function(file, ...) {
  message(paste("Reading",file))
  if ( tools::file_ext(file) %in% c('xlsx','xls') ) {
    data = readxl::read_xlsx(file,1,...)
  } else {
    data = readr::read_tsv(file, col_names=T, trim_ws=T, ...)  
  }
  return(data)
}

# Loads a data file with per-sample information on pathology, p53 IHC, etc
.setUpRxTablePerEndo<-function(predDf, func=c('max','mean'), verbose=T) {
  func = match.arg(func)
  if (verbose)
    message(paste0("Evaluating the ",func," risk per endoscopy"))
  
  if (func == 'max') {
    predDf = predDf %>% group_by(Endoscopy) %>% dplyr::mutate('n.samples'=length(Endoscopy),'Samples'=paste(Sample,collapse=',')) %>% 
      dplyr::summarise_if(is.numeric, list(max)) %>% dplyr::select(-matches('^Sample$'))
  } else {
    predDf = predDf %>% group_by(Endoscopy) %>% dplyr::mutate('n.samples'=length(Endoscopy), 'Samples'=paste(Sample,collapse=',')) %>%  
        dplyr::summarise_if(is.numeric, list(mean))
  }

  return(predDf)
}
