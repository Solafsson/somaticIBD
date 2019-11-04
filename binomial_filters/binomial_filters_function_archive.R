
## The functions below are for binomial filtering of somatic calls to remove germline variants and similar work.
## Almost everything comes from Tim Coorens.

#-------------------------------------------------
# Filters for various projects:
# Exact binomial test - germline variants
# Beta-binomial overdispersion - artefacts
# Binomial mixture model - noise/contamination/subclonal variants
# Tim Coorens - January 2019
#-------------------------------------------------

options(stringsAsFactors = F)

if(file.exists("/nfs/users/nfs_s/so11/phd/doIexist")) {
  ## Script is running on the farm
  lib_location="/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm3/"
  .libPaths( c(lib_location, .libPaths()) )
  library(ggplot2, lib.loc=lib_location)
  #library(cowplot, lib.loc=lib_location)
  require(VGAM, lib.loc=lib_location)
} else {
  ## Running the script locally - Testing
  library(ggplot2)
  library(cowplot)
  require(VGAM)
}


#-------------------------------------------------
# Functions
#-------------------------------------------------

#1. Exact binomial test
## Split into X and Y chromosomes and autosomes.
## NV: Reads supporting variant
## NR: Total reads
## Keep variants that are overdispersed.
## VAF in an individual sample, not the global vaf.

## bionom_mix is a wrapper for the other three functions.

## Filters out recurrent artifacts. Present in multiple/all samples at low vafs.

## rho cut-off of 0.1
## could estimate best value by doing a histogram of the log(rho) across all samples.

exact.binomial=function(gender,NV,NR,cutoff=-5){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting
  # variant and NR the matrix with total depth (samples as columns, mutations rows,
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector,
  # TRUE if mutation is likely to be germline.

  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal

  if(gender=="female"){
    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))

    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  qval = p.adjust(pval,method="BH")
  germline = log10(qval)>cutoff
  return(germline)
}


#2. Beta-binomial overdispersion filter


estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV,cutoff=0.1, binom.pval=F,pval.cutoff=0.05){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input.
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE

  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
    if(binom.pval){
      mu = sum(as.numeric(NV[k,]))/sum(as.numeric(NR[k,]))
      LL0 = sum(dbinom(as.numeric(NV[k,]),as.numeric(NR[k,]),prob=mu))
      LL1 = sum(dbetabinom(as.numeric(NV[k,]),as.numeric(NR[k,]),prob=mu,rho=rho_est[k]))
      pval[k] = (1-pchisq(2*(LL1-LL0),df=1)) / 2
    }
    if (k%%1000==0){
      print(k)
    }
  }
  if(binom.pval){
    qval=p.adjust(pval,method="BH")
    flt_rho=qval<=pval.cutoff&rho_est>cutoff
  }else{
    flt_rho=rho_est>=cutoff
  }
  return(flt_rho)
}

#3. Binomial Mixture Model

## Expectation step
estep = function(x,size,p.vector,prop.vector,ncomp){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood

  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }

  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

## Maximisation step
mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)

  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)
}

## EM algorithm
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities

  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]

  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]

    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}

    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]

  }
  if(!flag) warning("Didn’t converge\n")

  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC)
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()

  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers

    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}



## Some functions that I added:
## vafs is a numeric vector of vafs. All mutations with VAF<lowest will be excluded. Name is for the plot.
## The function returns a plot object. Needs to be printed.
plot_vaf_distribution <- function(vafs, lowest=0, name="Crypt") {
  thisCrypt <- vafs[vafs>lowest]
  median <- median(thisCrypt)
  p <- ggplot(data.frame(VAF=thisCrypt), aes(x=VAF)) + geom_histogram() + ggtitle(paste(name, "VAF distribution \n Median:", median)) + theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

## To stick in apply()
find_median_vaf <- function(vafs, lowest=0) {
  thisCrypt <- vafs[vafs>lowest]
  return(median(thisCrypt, na.rm=T))
}

Mutcount <- function(vaf_vect, lower=0) {
  count <- length(vaf_vect[vaf_vect>lower & !is.na(vaf_vect)])
  return(count)
}



