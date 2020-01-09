
To generate the report provided with this package you need to install pandoc version 1.12.3. If you do not also have LaTEX installed, follow the instructions from the pandoc link to install as well.

http://pandoc.org/installing.html


# Example Usage for Predicting Barrett's Progression

Using your own data and the qdnaseq.R script found in the example/ directory first process one or more BAM files from a single patient and output the files into a unique directory. Then run `predictRisk(path='.')`

```

library(BarrettsProgressionRisk)

# bamPath needs to contain one or more bam files for a patient
runQDNAseq(bamPath='.', outputPath='qdnaseq_output/')

segObj = segmentRawData(loadSampleInformation('example/endoscopy.xlsx'), 
  raw.data = 'example/qdnaseq.output.binSize50.rawReadCounts.txt', 
  fit.data = 'example/qdnaseq.output.binSize50.fittedReadCounts.txt', verbose=F)

pr = predictRiskFromSegments(segObj, verbose=F)

## Results

# get QC information
sampleResiduals(pr)

# plot raw data
plotSegmentData(pr)

# risk per sample
predictions(pr,'sample')

# risk per endoscopy
predictions(pr,'endoscopy')


# output the absolute risk CI per sample
absoluteRiskCI(pr, 'sample')

# output the absolute risk CI per endoscopy
absoluteRiskCI(pr, 'endoscopy')


# Get recommendations per endoscopy, given as either sequential integers or dates in the sample information loaded initially.

rx(pr)

```

# Report Generator

A script is included in inst/bin/generate_report.R. To write your own three parameters are needed:

1. The output directory where the raw and fitted files were output from QDNAseq. These must have been generated with a bin size of 50kb (default if the provided function is used).
2. A file containing per-sample p53 IHC and pathology information (see example 'endoscopy.xml'). This is optional.
3. The output directory for the html report.

```
library(rmarkdown)
library(knitr)

# bamPath needs to contain one or more bam files for a patient, the default binsize is 50. 
# This should not be changed without extensive testing of the data unless you are retraining the underlying model!  
BarrettsProgressionRisk::runQDNAseq(bamPath='.', outputPath=<path to qdnaseq output>,  binsize=50)


# qdnaseq.path=<path to qdnaseq output>
qdnaseq.path='examples/'
# info.file=<path to per sample p53 IHC/pathology file>
info.file = 'example/endoscopy.xlsx'

output.dir='~/tmp'

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
rmarkdown::render(rmd, params=list(path=path.expand(qdnaseq.path), info.file=info.file), 
                  output_dir=path.expand(output.dir), output_format='html_document',
                  intermediates_dir=path.expand(output.dir))

```


