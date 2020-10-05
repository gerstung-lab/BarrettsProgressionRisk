require(tidyverse)

be.model.fit<-function(model, s, kb, tile.size, 
                       tile.mean, arms.mean, tile.sd, arms.sd, 
                       cx.mean, cx.sd, per.pt.nzcoefs, cvRR, df = NULL, pconf = NULL) {

  if (is.null(pconf)) pconf = be_model$pred.confidence
  
  be.model <- list(
    fit = model, lambda = s, qdnaseq_kb = kb, tile.size = tile.size,
    tile.mean = tile.mean, arms.mean = arms.mean, tile.sd = tile.sd,
    arms.sd = arms.sd, cx.mean = cx.mean,  cx.sd = cx.sd, 
    nzcoefs = per.pt.nzcoefs, cvRR = cvRR, fit.data = df,
    pred.confidence = pconf
  )
  

  class(be.model) <- c('BEModel', class(be.model))
    
  return(be.model)
}
 



pi.hat<-function(x) exp(x)/(1+exp(x))


#' Times per sample are determined by the order of the sample as given by the demoFile, or by the order of the samples in the dataset.
#' @name rxRules
#' @return A table describing recommendation rules
#' @author skillcoyne
#' @export
rxRules<-
  tibble(
    'Rule' = c(1:4,'None'),
    'Rx' = c('Immediate RFA', 'Recheck in 6-12 months',
             'Recheck in 12-24 months','Regular surveillance in 3-5 years','N/A'),
    'Description' = c(
      'HGD or IMC diagnosis or more than one consecutive high risk predictions.',
      'One high risk prediction or an aberrant p53 IHC.',
      'One or more moderate risk predictions.',
      'Two or more consecutive low risk predictions.','')
  )



tileSamples<-function(obj, be.model=NULL, scale=T, MARGIN=2, verbose=T) {
  if (verbose) message(paste0('Scale tiled data: ',scale))
  
  if (is.null(be.model)) {
    be.model = BarrettsProgressionRisk:::be_model
  } else {
    if (length(be.model$arms.mean) != length(be.model$arms.sd)) stop('Arm means/sd do not match in length') 
    if (length(be.model$tile.mean) != length(be.model$tile.sd)) stop('Tile means/sd do not match in length') 
  }

  # Tile, scale, then merge segmented values into 5Mb and arm-length windows across the genome.
  segtiles = tileSegments(obj, size = be.model$tile.size, verbose=verbose)
  armtiles = tileSegments(obj, size='arms',verbose=verbose)
  if (scale & MARGIN == 2) {
    if (verbose) message('Scaling and centering per bin')
    for (i in 1:ncol(segtiles$tiles))
      segtiles$tiles[,i] = unit.var(segtiles$tiles[,i], be.model$tile.mean[i], be.model$tile.sd[i])
    for (i in 1:ncol(armtiles$tiles))
      armtiles$tiles[,i] = unit.var(armtiles$tiles[,i], be.model$arms.mean[i], be.model$arms.sd[i])
  } else if (scale & MARGIN == 1) {
    if (verbose) message('Scaling and centering per sample')
    tl = t(apply(segtiles$tiles, 1, scale, center=T, scale=T))
    colnames(tl) = colnames(segtiles$tiles)
    segtiles$tiles = tl
    
    tl = t(apply(armtiles$tiles, 1, scale, center=T, scale=T))
    colnames(tl) = colnames(armtiles$tiles)
    armtiles$tiles = tl
  }
    
  cx.score = scoreCX(segtiles$tiles,1)
  if (scale & MARGIN == 2) {
    cx.score = unit.var(cx.score, be.model$cx.mean, be.model$cx.sd)
  } else if (scale & MARGIN == 1) {
    cx.score = cx.score/sqrt(mean(be.model$cx.mean^2))
  }
  
  mergedDf = subtractArms(segtiles$tiles, armtiles$tiles)
  mergedDf = cbind(mergedDf, 'cx'=cx.score)
  
  retlist = list('tiles'=mergedDf, 'residuals'=NULL, 'per.sample.error'=NULL)
  
  if (!is.null(segtiles$error)) {
    pse = per.sample.error(mergedDf, segtiles$error, armtiles$error, be.model)
    retlist = list('tiles'=mergedDf, 'residuals'=pse$Xerr, 'per.sample.error'=pse$error)
  }
  
  return(retlist)  
}

per.sample.error<-function(df, seg.tile.error, arm.tile.error, be.model) {
  # get the bootstrap errors for coefficients  
  coef.error = .bootstrap.coef.stderr(be.model)
  
  # Errors are the per window, weighted mean error of all segments in the bin
  Xerr = cbind(seg.tile.error, arm.tile.error)
  Xerr_diag = apply(Xerr, 1, function(x) { diag(as.matrix(x^2)) })
  
  # covariance
  covB = cov(coef.error[,'jack.se',drop=F])[1]
  X = t(df[,coef.error$coef,drop=F])
  B = coef.error[,'1',drop=F] %>% data.frame
  
  Var_rr = apply(X,2,function(x) sum(x*covB*x)) + sapply(Xerr_diag, function(x) sum(B*x*B))
  perSampleError = sqrt( Var_rr )
  return(list('error'=perSampleError, 'residuals'=Xerr))
}


predictRisk<-function(merged.tiles, be.model=NULL, obj=NULL, info=NULL, verbose=T) {
  if (is.null(be.model)) {
    be.model = BarrettsProgressionRisk:::be_model
  } 
  
  if (is.null(obj) & is.null(info))
    stop("Need either a SegmentedSWGS or SampleInformation object")

  sample.info = info
  if (!is.null(obj)) {
    sample.info = obj$sample.info %>% dplyr::filter(Sample %in% (sampleResiduals(obj) %>% dplyr::filter(Pass) %>% dplyr::select(dplyr::matches('sample')) %>% pull))
  }
  
  sparsed_test_data = Matrix(data=0, nrow=nrow(merged.tiles$tiles),  ncol=ncol(merged.tiles$tiles),
                             dimnames=list(rownames(merged.tiles$tiles),colnames(merged.tiles$tiles)), sparse=T)
  for(i in colnames(merged.tiles$tiles)) sparsed_test_data[,i] = merged.tiles$tiles[,i]
  
  if (!is.null(merged.tiles$per.sample.error)) {
    perSampleError = tibble('Sample'=names(merged.tiles$per.sample.error),'Error'=merged.tiles$per.sample.error)
    perSampleError = left_join(perSampleError, sample.info, by='Sample') %>% dplyr::select('Sample','Error','Endoscopy')
  } else {
    perSampleError = NULL
  }
  
  # Predict and generate absolute probabilities
  RR = predict(be.model$fit, newx=sparsed_test_data, s=be.model$lambda, type='link')
  probs = BarrettsProgressionRisk:::pi.hat(RR)
  
  per.sample.preds = full_join(tibble('Sample'=rownames(probs), 'Probability'=round(probs[,1],2), 
                                      'Relative Risk'=RR[,1], 'Risk'=sapply(probs[,1], BarrettsProgressionRisk:::.risk, be.model$pred.confidence)), 
                              sample.info, by='Sample')
  
  per.endo.preds = BarrettsProgressionRisk:::.setUpRxTablePerEndo(per.sample.preds, 'max', verbose) %>% rowwise() %>% 
    mutate( Risk = BarrettsProgressionRisk:::.risk(Probability,be.model$pred.confidence))
  
  if (!is.null(merged.tiles$per.sample.error)) {
    perEndoError = BarrettsProgressionRisk:::.setUpRxTablePerEndo(perSampleError, 'max', F) %>% dplyr::select('Endoscopy','Error')
  } else {
    perEndoError = NULL
  }
  
  psp = list('per.endo'=per.endo.preds, 'per.sample'=per.sample.preds, 'segmented'=obj, 'per.sample.error'=perSampleError, 'per.endo.error'=perEndoError, 'tiles'=merged.tiles$tiles, 'tiles.resid'=merged.tiles$residuals, 'be.model'=be.model, 'sample.info'=sample.info)
  
  class(psp) <- c('BarrettsRiskRx', class(psp))
  return(psp)
}


#' Get the per-sample residuals calculated from the segmentation phase.
#' @name predictRiskFromSegments
#' @param SegmentedSWGS object
#' @return Predict segmented data using included model
#'
#' @author skillcoyne
#' @export
predictRiskFromSegments<-function(obj, be.model = NULL, verbose=T) {
  if (class(obj)[1] != 'SegmentedSWGS')
    stop("SegmentedSWGS object missing")
  
    if (nrow(sampleResiduals(obj)) == nrow(sampleResiduals(obj) %>% dplyr::filter(!Pass))) {
      warning('No samples passed QC, no predictions can be made.')
      return(NULL)
    }

  if (is.null(be.model)) {
    be.model = BarrettsProgressionRisk:::be_model
    message('Using internal glmnet model.')
  } else {
    warning("Using EXTERNAL glmnet model. Validation not provided.")
  }
    
  # Tile, scale, then merge segmented values into 5Mb and arm-length windows across the genome.
  binnedSamples = tryCatch({
    BarrettsProgressionRisk:::tileSamples(obj, be.model, scale=T, MARGIN=2, verbose=verbose)
  }, error = function(e) {
    msg = paste("ERROR tiling segmented data:", e)
    stop(msg)
  })
  
  # Predict  
  psp = predictRisk(binnedSamples, be.model, obj, verbose)
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
             'Risk.low'=sapply(low, .risk, psp$be.model$pred.confidence),
             'CI.high'=high,
             'Risk.high'=sapply(high, .risk, psp$be.model$pred.confidence))
  
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

  # TODO calculate this from the model!
  total = cbind.data.frame('NP'=43,'P'=45)
  cases = total[['P']]
  
  # so...not sure why I'm multiplying by 100 this way.  This is the same as saying cases/0.01
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
  adjustedRiskRx$Risk = sapply(adjustedRiskRx$Probability, .risk, be.model$pred.confidence)

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

  rules = .apply.rules(preds,riskBy, brr$be.model$pred.confidence)

  return(rules)
}

.apply.rules<-function(preds,riskBy,pred.confidence,time=c('date','numeric')) {
  time = match.arg(time)
  
  if (riskBy == 'Sample') warning('Rx rules are intended to be applied on a per-endoscopy basis, not per-sample.')

  time = 'date'  
  if (!class(preds[[riskBy]]) %in% 'Date' || is.numeric(preds[[riskBy]])) {
    time = 'numeric'
    preds[[riskBy]] = factor(preds[[riskBy]])
  }
  
  preds = preds %>% rowwise() %>%
    dplyr::mutate(Risk = .risk(Probability, pred.confidence)) %>%
    mutate(Risk = factor(Risk, levels=names(riskColors()), ordered=T))
  
  p53Col = grep('p53', colnames(preds), value=T, ignore.case=T)
  pathCol = grep('pathology', colnames(preds), value=T, ignore.case=T)
  
  # Consecutive after sorting by the selected column
  preds = preds %>% arrange(preds[[riskBy]])
  
  rules = tibble(`Time 1` = as.Date(x = integer(0), origin = "1970-01-01"), 
                 `Time 2` = as.Date(x = integer(0), origin = "1970-01-01"), 
                 Rule = integer())
  if (time == 'numeric') 
    rules = tibble(`Time 1` = double(), `Time 2` = double(), Rule = integer())
  
  if (riskBy == 'Sample')
    rules = add_column(rules,`Sample 1` = character(), `Sample 2` = character())  

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
    
    t1 = preds[[riskBy]][i] 
    t2 = preds[[riskBy]][(i+1)]
    if (time == 'date') {
      t1 = as.Date(t1); t2 = as.Date(t2)
    } else {
      t1 = as.integer(t1); t2 = as.integer(t2)
    }

    if (riskBy == 'Sample') {
      rules = add_row(rules, 'Time 1' = t1, 'Time 2' = t2, 'Rule' = rule, 'Sample 1' = preds[[riskBy]][i], 'Sample 2' = preds[[riskBy]][i+1]  ) 
    } else {
      rules = add_row(rules, 'Time 1' = t1, 'Time 2' = t2, 'Rule' = rule) 
    }

    if (i == 1 && i == nrow(preds)) break;
    
    # If the last pair resulted in an increased surveillance rx than this one should do as well?
    if ( i == nrow(preds) && rules[i,'Rule'] > rules[(i-1), 'Rule'] ) rules[i,'Rule'] = rules[(i-1), 'Rule'] 
    
    if (i == nrow(preds)) break;
  }
  rules = rules %>% dplyr::mutate(Rx = purrr::map_chr(Rule, .rule.rx))
  

  return(rules)
}

# Current rules
.rule.rx<-function(n) {
  return( subset(rxRules, Rule == n)$Rx )
}


.risk<-function(p, pred.confidence) {
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
      dplyr::summarise_if(is.numeric, list(max)) %>% dplyr::select(-dplyr::matches('^Sample$'))
  } else {
    predDf = predDf %>% group_by(Endoscopy) %>% dplyr::mutate('n.samples'=length(Endoscopy), 'Samples'=paste(Sample,collapse=',')) %>%  
        dplyr::summarise_if(is.numeric, list(mean))
  }

  return(predDf)
}
