
To generate the report provided with this package you need to install pandoc version 1.12.3. If you do not also have LaTEX installed, follow the instructions from the pandoc link to install as well.

http://pandoc.org/installing.html



# Example Usage

```
library(BarrettsProgressionRisk)

pr = predictRisk(path='.')

# Plot raw data
plotSegmentData(pr)

# output the risk table per sample
predictions(pr)

# Get recommendations per sample pair, assuming these were at different timepoint
rx(pr, 'demo_file.txt')

# Samples that failed QC
sampleNames(pr,F)
```

# Report Generator

A script is included in inst/bin/generate_report.R. To write your own three parameters are needed:

1. The output directory where the raw and fitted files were output from QDNAseq
2. A clinical file containing per-sample p53 IHC and pathology information (see example 'demo_file.txt'). This is optional.
3. The output directory for the PDF report.


```
library(rmarkdown)
library(knitr)

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
rmarkdown::render(rmd, params=list(path='.', info.file=NULL), output_dir='.')

```