args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2)
  stop("Usage: <directory containing QDNAseq files> <clinical info file>  <output dir> ")

suppressPackageStartupMessages( library(rmarkdown) )
suppressPackageStartupMessages( library(knitr) )

print(args)
print(length(args))

qdnaseq.path = args[1]
clin.file = args[2]
output.dir = args[3]

if (length(grep('raw', list.files(qdnaseq.path,pattern='txt'))) <= 0)
  stop(paste("No raw file found in",qdnaseq.path))

if(length(grep('corr|fitted', list.files(qdnaseq.path,pattern='txt'))) <= 0)
  stop(paste("No fitted file found in",qdnaseq.path))

if (!file.exists(clin.file)) 
  stop(paste("Clinical information file",clin.file,"does not exist"))


rawFiles = grep('raw',list.files(qdnaseq.path,'txt',full.names=T), value=T)
fittedFiles = grep('fitted',list.files(qdnaseq.path,'txt',full.names=T), value=T)

tmp.input = NULL
if (length(rawFiles) > 1) {
  if (length(rawFiles) != length(fittedFiles))
    warning(paste("Missing 1 or more raw or fitted files in",qdnaseq.path))
  
  raw.data = NULL; fitted.data = NULL
  for (file in rawFiles) {
    print(file)
    raw = data.table::fread(file)
    fit = data.table::fread(grep(sub('\\.raw.*', '', basename(file)), fittedFiles, value=T))
  
    if (ncol(raw) == 5) {
      colnames(raw)[5] = sub('\\.raw.*','',basename(file))
      colnames(fit)[5] = sub('\\.raw.*','',basename(file))
    }
    
    if (is.null(fitted.data)) {
      fitted.data = fit
      raw.data = raw
    } else {
      fitted.data = base::merge(fitted.data, fit, by=c('location','chrom','start','end'), all=T) 
      raw.data = base::merge(raw.data, raw, by=c('location','chrom','start','end'), all=T) 
    }
  }
  
  tmp.input = file.path(tempfile(tmpdir = output.dir, pattern='qdnaseq'), fsep = .Platform$file.sep)
  dir.create(tmp.input, recursive = T)

  write.table(raw.data[order(raw.data$chrom, raw.data$start),], row.names=F, quote=F, sep='\t', file=paste0(tmp.input, .Platform$file.sep, 'rawReadCounts.txt'))
  write.table(fitted.data[order(fitted.data$chrom, fitted.data$start),], row.names=F, quote=F, sep='\t', file=paste0(tmp.input,.Platform$file.sep,'fittedReadCounts.txt'))
  
  qdnaseq.path = tmp.input
}

qdnaseq.path = normalizePath(qdnaseq.path)
clin.file = normalizePath(clin.file)
output.dir = normalizePath(output.dir)

message(paste('Data path:', qdnaseq.path))
message(paste('Information file:', clin.file))
message(paste('Output directory:', output.dir))

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
message(rmd)

rmarkdown::render(rmd, params=list(path=qdnaseq.path, info.file=clin.file), 
                    output_dir=output.dir, output_format='html_document',
                  intermediates_dir=output.dir)

if (!is.null(tmp.input))
  unlink(tmp.input, recursive = T)  

message(paste("Report saved to: ", normalizePath(output.dir),'/RiskReport.html', sep=''))

