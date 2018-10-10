args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2)
  stop("Usage: <directory containing QDNAseq files> <clinical info file> <output dir OPT>")

qdnaseq.path = args[1]
clin.file = args[2]

output.dir = '.'
if (length(args == 3)) output.dir = args[3]

if (length(grep('raw', list.files(qdnaseq.path,pattern='txt'))) <= 0)
  stop(paste("No raw file found in",qdnaseq.path))

if(length(grep('corr|fitted', list.files(qdnaseq.path,pattern='txt'))) <= 0)
  stop(paste("No fitted file found in",qdnaseq.path))

if (!file.exists(clin.file)) {
  warning(paste("Clinical information file",clin.file,"does not exist"))
  clin.file = NULL
}

print(qdnaseq.path)
print(clin.file)
print(output.dir)


library(rmarkdown)
library(knitr)

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
rmarkdown::render(rmd, params=list(path=path.expand(qdnaseq.path), info.file=clin.file), output_dir=path.expand(output.dir))

message(paste("Report saved to: ", output.dir,'/RiskReport.pdf', sep=''))

