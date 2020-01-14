## THIS WOULD HAVE TO CHANGE IF THE INTERNAL MODEL IS CHANGED
test_that('internal data', {
  
  expect_equal(length(BarrettsProgressionRisk:::be_model), 14)
  
  expect_true('glmnet' %in% class(BarrettsProgressionRisk:::be_model$fit))
  
  expect_true(is.double(BarrettsProgressionRisk:::be_model$lambda))
  
  expect_equal(BarrettsProgressionRisk:::be_model$qdnaseq_kb, 50)
  
  expect_equal(BarrettsProgressionRisk:::be_model$tile.size, 5e6)
  
  expect_equal(dim(BarrettsProgressionRisk:::be_model$cvRR),c(75,3))
  
  expect_equal(dim(BarrettsProgressionRisk:::be_model$fit.data), c(773, 634))
  
})

test_that("loadSampleInfo", {
  
  si = tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c(1,1,2) )
  
  expect_error(loadSampleInformation(), label='No data')
  expect_warning(loadSampleInformation(si), label='Missing useful columns (not required).')
  
  si = tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','FOO','LGD') )

  expect(loadSampleInformation(si), ok=T)
})


test_that('prep raw data', {
  
  raw.data = readr::read_tsv( '../../example/qdnaseq.output.binSize50.rawReadCounts.txt', col_names=T, col_types = cols('chrom'=col_character()))
  fit.data = readr::read_tsv( '../../example/qdnaseq.output.binSize50.fittedReadCounts.txt', col_names=T, col_type = cols('chrom'=col_character()))
  
  blacklist = system.file('extdata', 'qDNAseq_blacklistedRegions.txt',package="BarrettsProgressionRisk")
  
  prepped = BarrettsProgressionRisk:::prepRawSWGS(raw.data,fit.data,readr::read_tsv(blacklist, col_types='cdd'),verbose=F)

  expect_true(is.list(prepped))
  expect_equal(names(prepped), c('data','sdevs','good.bins','window.depths.standardised','fit.data','cv.plot'))
  
  expect_equal(dim(prepped$data), c(51801,8))
  expect_equal(length(prepped$sdevs), 6)
  expect(length(prepped$good.bins) > 0, ok=T)
  expect_equal(dim(prepped$window.depths.standardised), c(51932,6))
  expect_equal(dim(prepped$fit.data), c(51932,10))
  
  expect_equal(colnames(prepped$data), c('chrom','start','A','B','C','D','E','F'))
  expect_equal(colnames(prepped$fit.data), c('location','chrom','start','end','A','B','C','D','E','F'))
  expect_equal(length(prepped$cv.plot), 6)
  expect_true(is.ggplot(prepped$cv.plot$A))
})     
    

test_that('segment data', {
  
  raw.data = readr::read_tsv( '../../example/qdnaseq.output.binSize50.rawReadCounts.txt', col_names=T, col_types = cols('chrom'=col_character())) %>% dplyr::select(-D,-E,-F)
  fit.data = readr::read_tsv( '../../example/qdnaseq.output.binSize50.fittedReadCounts.txt', col_names=T, col_type = cols('chrom'=col_character())) %>% dplyr::select(-D,-E,-F)
  
  expect_equal(dim(raw.data), dim(fit.data))
  
  si = loadSampleInformation(tibble( 'Sample'=c('Z','X','Y'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','NDBE','LGD') ))
  
  expect_error(segmentRawData(si, raw.data, fit.data, verbose=F), label='No matching samples in data file and sample information.')


  info = loadSampleInformation(tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','NDBE','LGD') ))

  
  
  segObj = segmentRawData(info, raw.data, fit.data, verbose=F)
  
  expect_true(class(segObj)[1] == 'SegmentedSWGS')

  expect_equal(sampleResiduals(segObj), segObj$residuals)
  expect_true(sampleResiduals(segObj)$Pass[1])

  expect_true(is.list(segObj$seg.plots))
  expect_true(is.ggplot(segObj$seg.plots$A))
  
  expect_equal(dim(segObj$seg.vals), c(304,8))
})
  

test_that('predict risk', {
  
  raw.data = readr::read_tsv( '../../example/qdnaseq.output.binSize50.rawReadCounts.txt', col_names=T, col_types = cols('chrom'=col_character())) %>% dplyr::select(-D,-E,-F)
  fit.data = readr::read_tsv( '../../example/qdnaseq.output.binSize50.fittedReadCounts.txt', col_names=T, col_type = cols('chrom'=col_character())) %>% dplyr::select(-D,-E,-F)

  info = loadSampleInformation(tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','NDBE','LGD') ))
  
  segObj = segmentRawData(info, raw.data, fit.data, verbose=F)
  prr = predictRiskFromSegments(segObj, verbose=F)
  
  expect_true(class(prr)[1] == 'BarrettsRiskRx')
  
  expect_equal(names(prr), c("per.endo","per.sample","segmented","per.sample.error","per.endo.error","tiles","tiles.resid","be.model"))
  
  expect_equal(dim(prr$per.endo.error), c(2,2))
  expect_equal(dim(prr$per.sample.error), c(3,3))
  
  preds = predictions(prr)
  
  expect_equal( length(which(preds$Probability > 0 & preds$Probability <= 1)), 3 )
  
  expect_equal( length(which(preds$Risk %in% names(riskColors()))), 3 )
  
  expect_equal( length(which(relativeRiskCI(prr,'sample')$Error > 0 & relativeRiskCI(prr,'sample')$Error <= 1)), 3 )
  
  expect_equal( grep('CI.RR.high', colnames(relativeRiskCI(prr, 'sample'))), 11)
  expect_equal( grep('CI.RR.low', colnames(relativeRiskCI(prr, 'sample'))), 10)
  
  ci = filter(absoluteRiskCI(prr, 'sample'), Sample == 'B')
  
  expect_true(ci$CI.low < ci$Probability)
  expect_true(ci$CI.high > ci$Probability)

  expect_equal(nrow(rx(prr)),2)
  expect_equal(rx(prr)$Rule[1], rx(prr)$Rule[2])
    
})
       
          
  