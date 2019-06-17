


test_that("loadSampleInfo", {
  
  si = tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c(1,1,2) )
  
  expect_error(loadSampleInformation())
  expect_warning(loadSampleInformation(si))
  
  si = tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','FOO','LGD') )

  expect(loadSampleInformation(si), ok = T)
})


test_that('prep raw data', {
  
  raw.data = readr::read_tsv( '../../example/qdnaseq.output.binSize15.rawReadCounts.txt', col_names=T, col_types = cols('chrom'=col_character()))
  fit.data = readr::read_tsv( '../../example/qdnaseq.output.binSize15.fittedReadCounts.txt', col_names=T, col_type = cols('chrom'=col_character()))
  
  si = loadSampleInformation(tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','NDBE','LGD') ))
  expect_error(segmentRawData(si, raw.data, fit.data), label='No matching samples in data file and sample information.')

  si = si %>% mutate(Sample = paste0('Sample.',Sample))

  prepped = BarrettsProgressionRisk:::.prepRawSWGS(raw.data,fit.data,readr::read_tsv('inst/extdata/qDNAseq_blacklistedRegions.txt', col_types='cdd'),F)

  expect_true(is.list(prepped))
  expect_equal(names(prepped), c('data','sdevs','good.bins','window.depths.standardised','fit.data','cv.plot'))
  
  expect_equal(dim(prepped$data), c(170574,5))
  expect_equal(length(prepped$sdevs), 3)
  expect(length(prepped$good.bins) > 0, ok=T)
  
  expect_equal(dim(prepped$window.depths.standardised), c(170681,3))
  
  expect_equal(dim(prepped$fit.data), c(170681,7))
  
  expect_equal(colnames(prepped$data), c('chrom','start','Sample.A','Sample.B','Sample.C'))
  expect_equal(colnames(prepped$fit.data), c('location','chrom','start','end','Sample.A','Sample.B','Sample.C'))
  
  expect_equal(length(prepped$cv.plot), 3)
  
  expect_true(is.ggplot(prepped$cv.plot$Sample.A))

})     
    
    
test_that('segment data', {
  raw.data = readr::read_tsv( '../../example/qdnaseq.output.binSize15.rawReadCounts.txt', col_names=T, col_types = cols('chrom'=col_character()))
  fit.data = readr::read_tsv( '../../example/qdnaseq.output.binSize15.fittedReadCounts.txt', col_names=T, col_type = cols('chrom'=col_character()))
  
  si = loadSampleInformation(tibble( 'Sample'=c('A','B','C'), 'Endoscopy'=c("2001/01/01","2001/01/01","2001/10/21"), 'GEJ.Distance' = c(1,1,3), 'P53 IHC' = c(0,0,0), 'Pathology'=c('NDBE','NDBE','LGD') ))
  
  expect_error(segmentRawData(si, raw.data, fit.data), label='No matching samples in data file and sample information.')
  
  info = si %>% mutate(Sample = paste0('Sample.',Sample)) %>% filter(Sample == 'Sample.A')

  
  segObj = segmentRawData(info, raw.data, fit.data,verbose = T) 
  
  expect_true(class(segObj)[1] == 'SegmentedSWGS')

  expect_equal(sampleResiduals(segObj), segObj$residuals)
  expect_true(sampleResiduals(segObj)$Pass)

  expect_true(is.list(segObj$seg.plots))
  expect_true(is.ggplot(segObj$seg.plots$Sample.A))
  
  expect_equal(dim(segObj$seg.vals), c(518,6))
  
  
  prr = predictRiskFromSegments(segObj, verbose = T)
  
  expect_true(class(prr)[1] == 'BarrettsRiskRx')
  
  expect_equal(dim(prr$per.endo.error), c(1,2))
  expect_equal(dim(prr$per.sample.error), c(1,3))
  
  expect_equal(predictions(prr)$Probability, 0.73)
  expect_equal(predictions(prr)$Risk, 'High')
  
  expect_equal( grep('CI.RR.high', colnames(relativeRiskCI(prr, 'sample'))), 11)
  expect_equal( grep('CI.RR.low', colnames(relativeRiskCI(prr, 'sample'))), 10)
  
  expect_equal(nrow(rx(prr)),1)
  expect_equal(rx(prr)$Rule, 2)
    
})
       
          
  