args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2)
  stop("Usage: <directory containing QDNAseq files> <output dir> <clinical info file OPT> ")

suppressPackageStartupMessages( library(rmarkdown) )
suppressPackageStartupMessages( library(knitr) )

print(args)
print(length(args))

qdnaseq.path = args[1]
output.dir = args[2]

clin.file = NULL
if (length(args) == 3) clin.file = args[3]

if (length(grep('raw', list.files(qdnaseq.path,pattern='txt'))) <= 0)
  stop(paste("No raw file found in",qdnaseq.path))

if(length(grep('corr|fitted', list.files(qdnaseq.path,pattern='txt'))) <= 0)
  stop(paste("No fitted file found in",qdnaseq.path))

if (!is.null(clin.file) && !file.exists(clin.file)) 
  warning(paste("Clinical information file",clin.file,"does not exist"))


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

options(warn = -1)
rmd = system.file('rmd','RiskReport.Rmd',package="BarrettsProgressionRisk")
rmarkdown::render(rmd, params=list(path=path.expand(qdnaseq.path), info.file=clin.file), 
                  output_dir=path.expand(output.dir), output_format='html_document',
                  intermediates_dir=path.expand(output.dir))

if (!is.null(tmp.input))
  unlink(tmp.input, recursive = T)  

message(paste("Report saved to: ", output.dir,'/RiskReport.html', sep=''))

