.titleCase<-function(x) {
  gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",x, perl=TRUE)
}


cvRR<-function(df, coefs) {
  apply( exp(t(df[, rownames(coefs)])*coefs[,1]), 1, function(x) {
    sd(x)/mean(x)
  })
}

#TBD
model.pred.confidence<-function(df) {
  
}



# Gets a temp cache directory, not really used yet
getcachedir<-function() {
  tm <- Sys.getenv(c('TMPDIR', 'TMP', 'TEMP'))
  d <- which(file.info(tm)$isdir & file.access(tm, 2) == 0)
  if (length(d) > 0)
    tm[[d[1]]]
  else if (.Platform$OS.type == 'windows')
    Sys.getenv('R_USER')
  else
    '/tmp'
}

# Renames 'get.chr.lengths'
chrInfo<-function(chrs =  c(1:22, 'X','Y'), prefix='chr', build='hg19', file=NULL) {
  local_file =  paste(build,'_info.txt',sep='')
  
  local_file = system.file("extdata", local_file, package="BarrettsProgressionRisk")
  
  if (is.null(file) & file.exists(local_file))
    file = local_file

  if (!is.null(file) && file.exists(file)) {
    #message(paste("Reading chromosome information for build ",build," from local file.", sep=''))
    chr.lengths = read_tsv(file, col_types = cols( chrom = col_character(),
                                                   chr.length = col_double(),
                                                   chr.cent = col_double(),
                                                   cent.gap = col_double(),
                                                   genome.length = col_double()))
  } else {
    chr.lengths = read_tsv(paste('http://genome.ucsc.edu/goldenpath/help/', build, '.chrom.sizes',sep=''), col_names = F, col_types = 'cd') %>% set_names(c('chrom','chr.length'))
    chr.lengths = chr.lengths %>% dplyr::filter(chrom %in% paste(prefix,chrs, sep='')) %>% 
      mutate(chrom = factor(chrom, levels=paste(prefix,chrs,sep=''))) %>%
      arrange(chrom)
    
    cytoband.url = paste('http://hgdownload.cse.ucsc.edu/goldenPath',build,'database/cytoBand.txt.gz',sep='/')
    cytoband.file = paste(.Platform$file.sep,'tmp', paste(build,basename(cytoband.url),sep='_'), sep='/')

    tryCatch({
      download.file(cytoband.url, cytoband.file, cacheOK = T)
    }, error = function(e)
      stop(paste("Could not download", cytoband.url, "\n", e))
    )
    cytobands = read_tsv(cytoband.file, col_names = F, col_types = 'cddcc') %>%
      set_names( c('chrom','start','end','band','attr') ) %>%
      mutate(chrom = factor(chrom, levels=paste(prefix,chrs, sep='')))

    centromeres = cytobands %>% dplyr::filter(attr == 'acen')

    chr.lengths = left_join(chr.lengths, centromeres %>% dplyr::group_by(chrom) %>% dplyr::summarise(
      chr.cent=mean(range(start, end)),
      cent.gap = (max(end)-min(start))/2
    ), by='chrom')
    
    chr.lengths$genome.length = cumsum(as.numeric(chr.lengths$chr.length))
  }

  if (is.null(file))
    file = paste(.Platform$file.sep, 'tmp', .Platform$file.sep, build, '_info.txt', sep='')

  write.table(chr.lengths, sep='\t', row.names=F, file=file)

  chr.lengths = chr.lengths %>% mutate(chr = sub(prefix, '', chrom)) %>% mutate(chr = factor(chr, levels=chrs, ordered=T))
    
  return(as_tibble(chr.lengths))
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

.bootstrap.coef.stderr<-function(be.model) {
  fitCoefs = as.data.frame(as.matrix(glmnet::coef.glmnet(be.model$fit, be.model$lambda)))[-1,,drop=F]
  fitCoefs = fitCoefs[which(fitCoefs != 0),,drop=F]

  be.model$nzcoefs = purrr::map(be.model$nzcoefs, function(x) {
    if (is_tibble(x)) return(x)
    x %>% as_tibble(rownames = 'coef')
  })
  
  cfs = unique(unlist(sapply(be.model$nzcoefs, function(x) x %>% dplyr::select(coef) )))
  loo.coefs = data.frame(matrix(ncol=0,nrow=length(cfs)))
  loo.coefs$coef = cfs

  for (pt in names(be.model$nzcoefs)) 
    loo.coefs = dplyr::full_join(loo.coefs, be.model$nzcoefs[[pt]], by='coef')
  
  colnames(loo.coefs)[-1] = c(1:length(be.model$nzcoefs))
  loo.coefs[is.na(loo.coefs)] = 0
  
  ch = tibble::as_tibble(fitCoefs, rownames='coef')
  loo.ch = dplyr::left_join(ch, loo.coefs, by='coef')

  jk<-function(x) {
    jk = bootstrap::jackknife(x,mean)
    jk$jack.se
  }
  
  ch$jack.se = apply(loo.ch[,-1], 1, jk)
  return(ch)
}


get.loc<-function(df) {
  locs = do.call(rbind.data.frame, lapply(colnames(df), function(x) unlist(strsplit( x, ':|-'))))
  colnames(locs) = c('chr','start','end')
  locs[c('start','end')] = lapply(locs[c('start','end')], function(x) as.numeric(as.character(x)))
  locs$chr = factor(locs$chr, levels=c(1:22,'X','Y'), ordered=T)
  as_tibble(locs)
}




non.zero.coef<-function(model=fitV, s=lambda) {
  cf = as.matrix(coef(model, s))
  cf[cf!=0,]
}

## Not used currently...
.mergeCountFiles<-function(path='.', raw.grep='raw.*read', corr.grep='corr|fitted', build='hg19', verbose=T) {
  chr.info = chrInfo(build=build)
  chr.info$chrom = sub('chr', '', chr.info$chrom)
  
  rawfiles = grep(raw.grep,list.files(path, 'txt'), value=T, ignore.case=T)
  
  pairlist = data.frame()
  for (name in sub('\\.raw.*','',rawfiles)) {
    pair = grep(name, list.files(path, 'txt', full.names=T), value=T)
    
    pairlist = rbind(pairlist, cbind('name'=name,
                                     'raw'=pair[which(grepl(raw.grep,pair,ignore.case=T))],
                                     'fitted'=pair[which(grepl(corr.grep,pair,ignore.case=T))]), stringsAsFactors=F)
  }
  
  merged.raw = NULL; merged.fitted = NULL
  for (i in 1:nrow(pairlist)) {
    name = pairlist[i,'name']
    dtR = data.table::fread(pairlist[i,'raw'])    
    colnames(dtR)[grep('loc|chr|start|end',colnames(dtR),invert=T)] = name
    dtF = data.table::fread(pairlist[i,'fitted'])    
    colnames(dtF)[grep('loc|chr|start|end',colnames(dtF),invert=T)] = name
    
    
    
    if (is.null(merged.raw)) {
      merged.raw = dtR
      merged.fitted = dtF
    } else {
      merged.raw = base::merge(merged.raw, dtR, by=c('location','chrom','start','end'), all=T) 
      merged.fitted = base::merge(merged.fitted, dtF, by=c('location','chrom','start','end'), all=T) 
    }
    
    if (verbose)
      message(paste( ncol(merged.raw)-4, ' samples merged. ', nrow(merged.raw), ' rows.', sep='' ))
    
  }
  
}
