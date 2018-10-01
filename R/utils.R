# Renames 'get.chr.lengths'
chrInfo<-function(chrs = paste('chr', c(1:22, 'X','Y'), sep=''), build='hg19', file=NULL) {
  local_file =  paste('inst/',build,'_info.txt',sep='')
  if (is.null(file) & file.exists(local_file))
    file = local_file

  if (!is.null(file) && file.exists(file)) {
    message(paste("Reading chromosome information for build ",build," from local file.", sep=''))
    chr.lengths = read.table(file, header=T, sep='\t', stringsAsFactors=F)
  } else {
    chr.lengths = read.table(paste('http://genome.ucsc.edu/goldenpath/help/', build, '.chrom.sizes',sep='') , sep='\t', header=F)
    colnames(chr.lengths) = c('chrom','chr.length')
    chr.lengths = subset(chr.lengths, chrom %in% chrs)
    chr.lengths$chrom = factor(chr.lengths$chrom, levels=chrs)
    chr.lengths = arrange(chr.lengths, chrom)

    cytoband.url = paste('http://hgdownload.cse.ucsc.edu/goldenPath',build,'database/cytoBand.txt.gz',sep='/')
    cytoband.file = paste('/tmp', paste(build,basename(cytoband.url),sep='_'), sep='/')

    tryCatch({
      download.file(cytoband.url, cytoband.file, cacheOK = T)
    }, error = function(e)
      stop(paste("Could not download", cytoband.url, "\n", e))
    )

    cytobands = read.table(cytoband.file, sep='\t', header=F)
    colnames(cytobands) = c('chrom','start','end','band','attr')
    cytobands$chrom = factor(cytobands$chrom, levels=chrs)

    centromeres = subset(cytobands, attr == 'acen')

    chr.lengths = cbind(chr.lengths,  centromeres %>% dplyr::group_by(chrom) %>% dplyr::summarise(
      chr.cent=mean(range(start, end)),
      cent.gap = (max(end)-min(start))/2
    ))
    chr.lengths$genome.length = cumsum(as.numeric(chr.lengths$chr.length))
  }

  # TODO -- not sure this will work installed on a Windows computer...
  if (is.null(file))
    file = paste('/tmp/', build, '_info.txt', sep='')

  write.table(chr.lengths, sep='\t', row.names=F, file=file)

  return(chr.lengths)
}


unit.var <- function(x, mean=NULL, sd=NULL, warn=F) {
  if ((is.null(mean) | is.null(sd)) || (is.na(mean) | is.na(sd))) {
    warning("No mean or sd provided.")
    if (length(x) == 1 | length(which(is.na(x))) == length(x) | sd(x, na.rm=T) == 0) {
      if (warn) warning("Unit normalization can't be performed with less than 2 samples or SD was 0")
      return(x)
    } else {
      return( (x-mean(x,na.rm=T))/sd(x,na.rm=T) )
    }
  } else {
    uv = (x-mean)/sd
    if (sd == 0) uv = 0
    return(uv)
  }
}



##Requires:
### medianFilter
getMad <- function(x,k=25){

  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]

  #Calculate runMedian - smoothing scatterplot
  runMedian <- medianFilter(x,k)

  dif <- x-runMedian # why subtract the smoothed values from the observations before taking the MAD?
  SD <- mad(dif) # Median Absolute Deviation

  return(SD)
}

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1

  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }

  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)
}
