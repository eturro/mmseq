readmmseq <- function( mmseq_files = grep("gene|identical", dir( pattern="\\.mmseq$" ), value=TRUE, invert=TRUE),
		 sample_names = sub("\\.mmseq","", mmseq_files), normalize=TRUE, uhfrac=NULL, partition="", common_set=NULL) {

	cat( "Using mmseq files:\n", paste("\t", mmseq_files, sep="", collapse="\n"), "\n") 
	cat( "Using sample names:\n", paste("\t", sample_names, sep="", collapse="\n"), "\n") 
	if( length(mmseq_files) == 0 || length(sample_names)  != length(mmseq_files) ) {
		cat( "ERROR:\n\tplease define mmseq files and sample names for each file.\n" )
		return(NULL)
	}
	
	# read in mmseq data
	mu = c()
	sd = c()
	mcse = c()
  iact = c()
	effective_lengths = c()
	true_lengths = c()
  n_reads <- c()
  unique_hits <- c()
  ntranscripts<- c()
  observed<- c()
  dup_sample_names <- c()
  percentiles <- list()
  percentiles_proportion <- list()
	cat( "Reading mmseq output files\n" )
	for( i in 1:length(mmseq_files) ) {
		dat = read.table( mmseq_files[i], header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t", quote="")
    if(! any(colnames(dat)=="log_mu") ||
       ! any(colnames(dat)=="sd") || 
       ! any(colnames(dat)=="mcse") || 
       ! any(colnames(dat)=="iact") || 
       ! any(colnames(dat)=="unique_hits") ||
       ! any(colnames(dat)=="ntranscripts") ||
       ! any(colnames(dat)=="observed") ||
       ! any(colnames(dat)=="effective_length") ||
       ! any(colnames(dat)=="true_length")) {
      stop("Colnames must contain log_mu, sd, mcse, iact, effective_length, true_length, ntranscripts, observed and unique_hits.")
    }
    spl = strsplit(partition, "|", fixed=TRUE)[[1]]
    if(length(spl)==0) spl=".*"
    if(length(spl)>0 && any(spl=="")) {
      stop("Cannot split using empty string.")
    }
    x=grep(spl[1], dat[,grepl("id", colnames(dat))]) 
    if(length(x)==0) spl=".*" # not a sample to partition

    for(p in spl) {
      cat("\t", mmseq_files[i], "(",p, ")","\t",sep="")
      dup_sample_names <- c(dup_sample_names, sample_names[i])
      if(p != ".*") dup_sample_names[length(dup_sample_names)] = paste(dup_sample_names[length(dup_sample_names)], p,sep="")

      x=grepl(p, dat[,grepl("id", colnames(dat))])
      if(!is.null(common_set)) {
        if(p != ".*") {
        #  x = x & sub(p, "", dat[,grepl("id", colnames(dat))]) %in% common_set
          x = match(common_set, sub(p, "", dat[, grepl("id", colnames(dat))]))
        } else {
        #  x = x & dat[,grepl("id", colnames(dat))] %in% common_set
          x = match(common_set, dat[,grepl("id", colnames(dat))])
        }
      }

      mu <- cbind(mu, dat$log_mu[x])
      sd <- cbind(sd, dat$sd[x])
      mcse <- cbind(mcse, dat$mcse[x])
      iact<- cbind(iact, dat$iact[x])
      unique_hits <- cbind(unique_hits, dat$unique_hits[x])
      ntranscripts<- cbind(ntranscripts, dat$ntranscripts[x])
      observed <- cbind(observed, dat$observed[x])
      effective_lengths <- cbind(effective_lengths, dat$effective_length[x])
      true_lengths <- cbind(true_lengths, dat$true_length[x])

      if(length(grep("percentiles[0-9]", names(dat))) == 1) {
        percentiles[[i]] <- do.call(rbind, lapply(dat[[grep("percentiles[0-9]", names(dat))]][x], function(y) {
          as.numeric(strsplit(y, ",")[[1]]) 
        }))
        colnames(percentiles[[i]])=strsplit(gsub("percentiles","", grep("percentiles[0-9]", names(dat), value=TRUE)), ",")[[1]]
      }

      if(length(grep("percentiles_proportion[0-9]", names(dat))) == 1) {
        percentiles_proportion[[i]] <- do.call(rbind, lapply(dat[[grep("percentiles_proportion[0-9]", names(dat))]][x], function(y) {
          as.numeric(strsplit(y, ",")[[1]]) 
        }))
        colnames(percentiles_proportion[[i]])=strsplit(gsub("percentiles_proportion","", grep("percentiles_proportion[0-9]", names(dat), value=TRUE)), ",")[[1]]
      }



      if(i==1) {
        rownames(mu) = rownames(sd) = rownames(mcse)= rownames(iact) = rownames(unique_hits)=rownames(ntranscripts)=rownames(observed)=rownames(true_lengths)=rownames(effective_lengths)=
          sub(partition, "", dat[x,grepl("id", colnames(dat))])
      } else {
        if(! all(rownames(mu) == sub(partition, "", dat[x,grepl("id", colnames(dat))]) )) {

          stop("Attempting to join incompatible mmseq files (different features).\n")
        }
      }
      ta <- readLines(mmseq_files[i], n=100) # header should be no more than 100 lines
      ta <- strsplit(ta[grep("# Mapped fragments", ta)], " ")[[1]]
      n_reads <- c(n_reads, as.numeric(ta[length(ta)]))
      cat( "found", n_reads[length(n_reads)], "aligned reads (or read pairs)\n" )
    }
  }
  colnames(mu)=colnames(sd)=colnames(mcse)=colnames(iact)=colnames(unique_hits)=colnames(ntranscripts)=colnames(observed)=colnames(effective_lengths)=colnames(true_lengths)=dup_sample_names
  names(n_reads) <- dup_sample_names

  if(!all(is.na(true_lengths)) && any(apply(true_lengths,1, function(x) { !all(x==x[1])}))) {
    warning("True transcript lengths differ across files.\n")
  }

  # un-log the estimates 
  cat( "Unlogging and unstandardising the estimates\n" )
  counts = exp(1)^(mu)
  
  # multiply by library size
  counts= sapply( colnames(counts), function(x) { counts[,x]*n_reads[x] })
  rownames(counts)=rownames(mu)

  # multiply by transcript length and divide by 10^9
  counts = sapply( colnames(mu), function(x) { counts[,x] = (counts[,x]*effective_lengths[,x])/(1e3*1e6) } )
  counts[is.na(counts)] = 0 # the ones that had no expression in MMSEQ
	
  if(normalize && ncol(mu)>1) {
    cat("Normalizing expression estimates with log factors:\n")
    if(is.null(uhfrac)) {
      minfrac = max(0.2, (ncol(mu) - floor(ncol(mu)^2/160))/ncol(mu))
    } else {
      minfrac=uhfrac
    }
    cat("Min frac: ", minfrac, "\n")
    uh = apply(unique_hits,1,function(x){sum(x>=1) >= minfrac*length(x) })
    sf <- apply(mu[uh,] - rowMeans(mu[uh,], na.rm=TRUE) ,2, function(x){median(x[is.finite(x)])})
    for(i in 1:length(sf)) {
      mu[,i] <- mu[,i] - sf[i]
      cat("\t", dup_sample_names[i],"\t", signif(sf[i],4),"\n")
    }
    cat("Normalizing counts with factors:\n")
    sf <- apply(log(counts[uh,]) - rowMeans(log(counts[uh,]), na.rm=TRUE), 2, function(x){median(x[is.finite(x)])})
    for(i in 1:length(sf)) {
      counts[,i] <- counts[,i]/exp(sf)[i]
      cat("\t", dup_sample_names[i],"\t", signif(exp(sf[i]),4),"\n")
    }
  }
  
  return(list(log_mu=mu, sd=sd, mcse=mcse, iact=iact, effective_length=effective_lengths, true_length=true_lengths, unique_hits=unique_hits, ntranscripts=ntranscripts[,1], observed=observed, counts=counts, percentiles=percentiles, percentiles_proportion=percentiles_proportion))
}

polyclass <- function(files, prior=NULL) {
  if(length(files)<2) {
    stop("Specify at least two mmdiff files")
  }
  n=length(files)+1
  if(is.null(prior)) {
    warning("Assuming flat prior across models")
    prior <- rep(1/n,n)
  } 
  if(length(prior) != n || !isTRUE(all.equal(sum(prior),1))) {
    stop("Specify a vector of ", n, " prior probabilities adding up to 1")
  }
  print("Note: make sure model 0 is the same in each comparison (only model 1 should differ)")

  res <- read.table(files[1],header=1)
  postm <- matrix(nrow=nrow(res), ncol=n)
  colnames(postm) <- paste("postprob_model",0:(n-1),sep="")
  bf <- matrix(nrow=nrow(res), ncol=n) # relative to model 1
  bf[,1] <- 1
  for(j in 1:length(files)) {
    bf[,j+1] <- read.table(files[j],header=1)$bayes_factor
  }

  for(i in 1:nrow(res)) {
    summ <- sum(bf[i,]*prior)
    postm[i,] <- bf[i,]*prior/summ
    if(sum(is.infinite(bf[i,]))>1) {
      warning(paste("More than one BF inf", i))
    } else if(sum(bf[1,]==0)>1) {
      warning(paste("More than one BF zero", i))
    } else if(sum(is.infinite(bf[i,]))==1) {
      postm[i,is.infinite(bf[i,])] <- 1
    }
  }
  w=grep("posterior_prob|bayes_factor|alpha|beta|eta", colnames(res))
  res <- res[,-w]
  res <- cbind(res,postm)
  return(res)
}


is.expressed <- function(hits_file, transcript_mmseq_file, unique_hits_thres=1) {

  if(system("hitstools", intern=FALSE, ignore.stderr=TRUE)>1) { 
    stop("Please make sure the hitstools tool is in your path")
  }

  gi=system(paste("hitstools header ", hits_file,  " | grep GeneIsoforms"), intern=TRUE)
  gi.l <- lapply(gi, function(x) { 
    sp=strsplit(x,"\t")[[1]]
    spl=sp[3:length(sp)]
    return(spl)
  })
  names(gi.l) <- sapply(gi, function(x) { strsplit(x,"\t")[[1]][2] })

  ig.l <- as.list(unlist(lapply(gi, function(x) { 
    sp=strsplit(x,"\t")[[1]]
    spl=rep(sp[2], times=length(sp)-2)
    names(spl)= sp[3:length(sp)]
    return(spl)
  })))

  tf=read.table(transcript_mmseq_file,header=1)
  log_mu_thres=min(tf$log_mu[tf$unique_hits>=unique_hits_thres])

  res.t=cbind(sapply(tf$feature_id,as.character),as.numeric(tf$log_mu>log_mu_thres))

  res=data.frame(gene=unlist(ig.l[match(res.t[,1],names(ig.l))]), transcript=res.t[,1], expressed=res.t[,2])

  res.spl=split(res,res$gene)

  res.g=lapply(res.spl, function(x) { 
    any(x$expressed==1)
  })
  
  return(list(transcript=cbind(sapply(res$transcript,as.character),sapply(res$expressed,as.character)), 
              gene=cbind(names(res.g), as.numeric(unlist(res.g)))))

}


ll.func <- function(pars, x=x,e=e, M=M, P=P) {
  N=length(x)
  C=ncol(P)
  V=ncol(M)
  par.alpha=pars[1]
  par.beta=pars[2:(V+1)]
  if(C==1) {
    par.gamma <- matrix(0,nrow=1,ncol=1)
    par.sigma2 <- matrix(pars[length(pars)],ncol=1)
  } else {
    par.gamma <- matrix(c(0,pars[(V+2):(V+C)]), ncol=1)
    par.sigma2 <- matrix(pars[(V+C+1):length(pars)],ncol=1)
  }
  Lambda = diag(1/as.vector(P %*% par.sigma2 + e))
  y=(x - rep(par.alpha,N) - M %*% par.beta - P %*% par.gamma) 
  return(-( rep(1,N) %*% diag(log(1/diag(Lambda))) %*% rep(1,N) + 
    t(y) %*% Lambda %*% (y)
  ))
}




mmseq.lrt <- function(x,e,M,P) {
  N=length(x)
  C=ncol(P)
  V=ncol(M)

  pars=c(1,rep(1,V), rep(1,C-1), rep(1,C))

  op = optim(pars, ll1.func, x=x, e=e, M=M, P=P, method="L-BFGS-B", control=list(fnscale= -1, trace=1), lower=c(rep(-Inf, length(pars)-2),0,0), upper=rep(Inf,length(pars)))

  mle1 <- list(value=exp(op$value),
    convergence=op$convergence,
    alpha=op$par[1],
    beta=op$par[2:(V+1)],
    gamma=c(0, op$par[(V+2):(V+C)]),
    sigma2=op$par[(V+C+1):(V+2*C)])

  P <- matrix(1, nrow=N,ncol=1)
  C=ncol(P)
  pars=c(1, rep(1,V), rep(1,C-1), 1)

  op = optim(pars, ll.func, x=x, e=e, M=M, P=P, method="L-BFGS-B", control=list(fnscale= -1, trace=1), lower=c(rep(-Inf, length(pars)-1),0), upper=rep(Inf,length(pars)))

  mle0 <- list(value=exp(op$value),
    convergence=op$convergence,
    alpha=op$par[1],
    beta=op$par[2:(V+1)],
    gamma=c(0),
    sigma2=op$par[(V+C+1):length(op$par)])

  return(list(stat=mle0$value/mle1$value, mle1=mle1, mle0=mle0))
}






ll0 <- function(x,e,M,P) {
  ll0.func <- function(pars, x=x, e=e, M=M) {
    N=length(x)
    V=ncol(M)
    par.alpha=pars[1]
    par.beta=pars[2:(V+1)]
    par.sigma2=pars[V+2]
    Lambda = diag(1/as.vector(P %*% par.sigma2 + e))
    y=(x - rep(par.alpha,N) - M %*% par.beta - P %*% par.gamma) 
    return(-.5 *( rep(1,N) %*% diag(log(1/diag(Lambda))) %*% rep(1,N) + 
      t(y) %*% Lambda %*% (y)
    ))
  }

  pars=c(rep(1,1),rep(1,ncol(M)), rep(1,ncol(P)-1), rep(1,ncol(P)))

  op = optim(pars, ll1.func, x=x, e=e, M=M, P=P, method="L-BFGS-B", control=list(fnscale= -1, trace=1), lower=c(rep(-Inf, length(pars)-2),0,0), upper=rep(Inf,length(pars)))

  N=length(x)
  C=ncol(P)
  V=ncol(M)
  opt <- list(value=op$value,
    convergence=op$convergence,
    alpha=op$par[1],
    beta=op$par[2:(V+1)],
    gamma=c(0, op$par[(V+2):(V+C)]),
    sigma2=op$par[(V+C+1):(V+2*C)])
  return(opt)
}








