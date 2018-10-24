
To generate the report provided with this package you need to install pandoc version 1.12.3. If you do not also have LaTEX installed, follow the instructions from the pandoc link to install as well.

http://pandoc.org/installing.html


# Example Usage for Predicting Barrett's Progression

Using your own data and the qdnaseq.R script found in the example/ directory first process one or more BAM files from a single patient and output the files into a unique directory. Then run `predictRisk(path='.')`

```

library(BarrettsProgressionRisk)

# bamPath needs to contain one or more bam files for a patient
runQDNAseq(bamPath='.', outputPath='qdnaseq_output/')

pr = predictRisk(path='qdnaseq_output/')

## Results

# get QC information
sampleResiduals(pr)

# plot raw data
plotSegmentData(pr)

# output the absolute risk per sample
predictions(pr)

# output the absolute risk CI per sample
absoluteRiskCI(pr)

# Get recommendations per sample pair (this assumes sequential timepoints per sample), if p53 IHC or pathology are available in a per-sample tab-separated file it can be included in this function call with rx(pr, myFile.txt)
rx(pr)


```

# Report Generator

A script is included in inst/bin/generate_report.R. To write your own three parameters are needed:

1. The output directory where the raw and fitted files were output from QDNAseq
2. A clinical file containing per-sample p53 IHC and pathology information (see example 'demo_file.txt'). This is optional.
3. The output directory for the html report.

```
library(rmarkdown)
library(knitr)

qdnaseq.path=<path to qdnaseq output>
info.file=<path to per sample p53 IHC/pathology file>

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
rmarkdown::render(rmd, params=list(path=path.expand(qdnaseq.path), info.file=clin.file), 
                  output_dir=path.expand(output.dir), output_format='html_document',
                  intermediates_dir=path.expand(output.dir))

```


