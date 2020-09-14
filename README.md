
To generate the report provided with this package you need to install pandoc version 1.12.3. If you do not also have LaTEX installed, follow the instructions from the pandoc link to install as well.

http://pandoc.org/installing.html


# Example Usage for Predicting Barrett's Progression 

## Using your own data

1. To process your own BAM files, ensure you have installed QDNAseq from Bioconductor. Recommend using the default binsize and minMapQ parameters as set in order to make your data comparable to the included model.
https://bioconductor.org/packages/release/bioc/html/QDNAseq.html

```
library(BarrettsProgressionRisk)

# Provide either a single bam file as an argument bam=<my file> or a path to one or more bam files path=<my bam path>

runQDNAseq(bam=<my file>, outputPath=<qdnaseq output path>)
```

2. Load the sample information for your samples.

```
# See the example of sample information table in example/endoscopy.xlsx 

info = loadSampleInformation(<my info table>)
```

3. Segment the fitted and raw QDNAseq files generated above.

```
segObj = segmentRawData(info, <outputPath/rawReadCounts.txt>, <outputPath/fittedReadCounts.txt>, verbose=F) 

```

4. Predict the risk for the segmented data.

```
sampleRisk = predictRiskFromSegments(segObj, verbose=F)

```

## Example using provided fitted and raw qdnaseq files:

```
library(BarrettsProgressionRisk)

# Load example data set
data(package='BarrettsProgressionRisk',ExampleQDNAseqData)

# bamPath needs to contain one or more bam files for a patient, the default binsize is 50. 
# This should not be changed without extensive testing of the data unless you are retraining the underlying model!  

# example of fitted data output from QDNAseq function call: runQDNAseq(bamPath='.', outputPath='<my path>/', binsize=50)
print(fit.data)
# example of raw data output from QDNAseq function call
print(raw.data)

# example of sample information table
print(info)

# Segment the fitted and raw data
segObj = segmentRawData(loadSampleInformation(info), raw.data, fit.data, verbose=F) 

# Predict risks from segmented data
pr = predictRiskFromSegments(segObj, verbose=F)

## -- Results  -- ##

# get QC information
sampleResiduals(pr)

# plot raw data
plotSegmentData(pr)

# Per sample mountain plots with annotations for coefficients
copyNumberMountainPlot(pr, annotate=T)

### Per sample classifications/probabilities

# risk per sample for samples that pass QC (see sampleResiduals(...))
predictions(pr,'sample')

# output the absolute risk CI per sample
absoluteRiskCI(pr, 'sample')

# Plot the sample classifications by Endoscopy and location 
patientRiskTilesPlot(pr)

### Per endoscopy classifications/probabilities

# risk per endoscopy
predictions(pr,'endoscopy')

# output the absolute risk CI per endoscopy
absoluteRiskCI(pr, 'endoscopy')

# Plot the endoscopy predictions over time, with risk classifications and confidence intervals based on the table output by absoluteRiskCI(...)
patientEndoscopyPlot(pr)

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
# BarrettsProgressionRisk::runQDNAseq(bamPath='.', outputPath=<path to qdnaseq output>,  binsize=50)


qdnaseq.path=system.file('extdata/example',package="BarrettsProgressionRisk")
info.file = system.file('extdata/example','endoscopy.xlsx',package="BarrettsProgressionRisk")
output.dir='~/tmp'

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
rmarkdown::render(rmd, params=list(path=path.expand(qdnaseq.path), info.file=info.file), 
                  output_dir=path.expand(output.dir), output_format='html_document',
                  intermediates_dir=path.expand(output.dir))

```


