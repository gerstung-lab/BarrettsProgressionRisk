.titleCase<-function(x) {
  gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",x, perl=TRUE)
}


cvRR<-function(df, coefs) {
  apply( exp(t(df[, rownames(coefs)])*coefs[,1]), 1, function(x) {
    sd(x)/mean(x)
  })
}


model.pred.confidence<-function(df) {
  ft.fun<-function(NP,P) {
    st.tb = df %>% dplyr::select(Patient, Status) %>% distinct %>% group_by(Status) %>% tally %>% spread(Status,n)
    f = chisq.test(rbind(cbind(NP,P),st.tb))
    cbind.data.frame('p.value'=round(f$p.value, 4))
  }
  
  qt = df %>% group_by(quants, Status) %>% dplyr::summarise(n=length(Status) ) %>% spread(Status, n) %>% ungroup %>%
    left_join( df %>% dplyr::group_by(quants) %>% dplyr::summarise ('mn'=mean(Probability), 'sd'=sd(Probability) ), by='quants')
  
  pred.confidence = qt %>% dplyr::group_by(quants) %>% 
    dplyr::mutate( 'perc'=P/sum(NP,P), 'p.value'=ft.fun(NP,P)$p.value, 'conf'=ifelse(p.value < 0.05, '*', '') ) %>% 
    separate(quants, c('r1','r2'), ',', remove = F) %>% 
    mutate_at(vars('r1','r2'), list(sub), pattern='\\[|\\]|\\(', replacement='') %>% mutate_at(vars(r1,r2), as.double) %>%  
    dplyr::mutate(Risk = ifelse(round(perc,1)<.5, 'Low','High'))
  
  pred.confidence = bind_cols(pred.confidence, 
                              data.frame(ci.low=qbeta(0.025, shape1=pred.confidence$P+.5, shape2 = pred.confidence$NP+.5),
                                         ci.high=qbeta(0.975, shape1=pred.confidence$P+.5, shape2 = pred.confidence$NP+.5)))
  
  pred.confidence %>% mutate(Risk = case_when( r2 == 0.5 | r1 == 0.3 ~ 'Moderate', TRUE ~ Risk ))
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
chrInfo<-function(chrs=c(1:22, 'X','Y'), prefix='chr', build='hg19', file=NULL) {
  local_file =  paste(build,'_info.txt',sep='')
  
  local_file = system.file("extdata", local_file, package="BarrettsProgressionRisk")
  tmp_file = paste(.Platform$file.sep, 'tmp', .Platform$file.sep, build, '_info.txt', sep='')
  
  if (is.null(file) & file.exists(local_file) & file.size(local_file) >= 1000) {
    file = local_file
  } else if (file.exists(tmp_file) & file.size(tmp_file) >= 1000) {
    file = tmp_file
  }

  if (!is.null(file) && file.exists(file)) {
    message(paste0("Reading chromosome information for build ",build," from ",file))
    
    chr.lengths = read.table(file, header = T, sep='\t', colClasses = c(character(), numeric(), numeric(), numeric(), numeric()), stringsAsFactors = F) %>% as_tibble()
    
  } else {
    chr.lengths = read.table(paste('http://genome.ucsc.edu/goldenpath/help/', build, '.chrom.sizes',sep=''), header = F, sep='\t', colClasses = c(character(), numeric()), stringsAsFactors = F) %>% 
      as_tibble() %>% set_names(c('chrom','chr.length'))
    
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
    cytobands = read.table(cytoband.file, header=F, colClasses = c(character(), numeric(), numeric(), character(), character()), stringsAsFactors = F ) %>% 
      as_tibble() %>% set_names( c('chrom','start','end','band','attr') ) %>%
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



generate.internal.be.model<-function(model.dir, saveObj=F) {
  select.alpha = '0.9'
  if (length(grep('loo_0.9|model_data|all.pt.alpha', list.files(model.dir))) < 3) stop(paste0("Missing required files in ", model.dir))
  
  load(paste0(model.dir,'/loo_0.9.Rdata'),verbose=F)
  rm(plots,performance.at.1se, fits, coefs)
  load(paste0(model.dir, '/model_data.Rdata'), verbose=F)
  load(paste0(model.dir, '/all.pt.alpha.Rdata'), verbose=F)
  load(paste0(model.dir, '/all.pt.alpha.Rdata'), verbose=F)

  orig.labels = tibble::enframe(labels) %>% mutate(SampleId = row_number())
  
  dysplasia.df = dysplasia.df[orig.labels$name,]
  rownames(dysplasia.df) = orig.labels$SampleId
  
  fitV = models[[select.alpha]]
  lambda = performance.at.1se[[select.alpha]]$lambda
  coefs = coefs[[select.alpha]]
  
  coef_cv_RR = tibble::enframe( BarrettsProgressionRisk:::non.zero.coef(fitV, lambda)[-1], name='label') %>% 
    dplyr::rename(coef = 'value') %>%
    mutate(cvRR = BarrettsProgressionRisk:::cvRR(dysplasia.df, coefs)[label]) %>% arrange(desc(cvRR))

  cuts = seq(0,1,0.1)  
  
  ids = pg.samp %>% dplyr::select(Hospital.Research.ID, Patient) %>% distinct %>% 
    dplyr::filter(Hospital.Research.ID %in% names(nzcoefs))

  nzcoefs = nzcoefs[ids$Hospital.Research.ID]
  names(nzcoefs) = as.character(ids$Patient)
  nzcoefs = purrr::map(nzcoefs, function(x) as_tibble(x, rownames = 'coef'))

  cxPredictions = pg.samp %>% 
    dplyr::select(Patient, Status, matches('Endoscopy'), Pathology, matches('Age|Sex|Gender'), Block, Samplename,  matches('Prediction|Probability'), RR) %>% 
    left_join(orig.labels, by=c('Samplename' = 'name') ) %>% dplyr::select(-Samplename, -value) %>% 
    dplyr::rename_at(vars(matches('Prediction')), list(~sub('Probability',.))) %>%
    dplyr::rename(Samplename = 'SampleId', `Relative Risk` = 'RR') %>% 
    dplyr::select(Patient, Samplename, Status, everything()) %>% 
    mutate(quants = cut(Probability, breaks=cuts, include.lowest = T))
  
  pred.conf = BarrettsProgressionRisk:::model.pred.confidence(cxPredictions)  
  
  
  be_model = BarrettsProgressionRisk:::be.model.fit(fitV, lambda, 5e6, z.mean, z.arms.mean, z.sd, z.arms.sd, mn.cx, sd.cx, nzcoefs, coef_cv_RR, pred.conf)  
  
  if (saveObj) saveRDS(be_model, file='R/sysdata.rda', version=2)
  
  return(be_model)
}
 


