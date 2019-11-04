if(file.exists("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt")) {
  lib_location="/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm3/"
  library(ape)
  library(parallel)
  .libPaths( c(lib_location, .libPaths()) )
  load_lib=function(){
  if(system("uname",intern = TRUE)=="Darwin"){
    dyn.load("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/phylogenics/treemut.Darwin.so")
  }else{
    dyn.load("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/phylogenics/treemut.so")
  }
}


} else {
  library(ape)
  library(parallel)

  load_lib=function(){
  if(system("uname",intern = TRUE)=="Darwin"){
    dyn.load("/Users/so11/phd/so11_nfs/somatic_ibd_p1/phylogenics/treemut.Darwin.so")
  }else{
    dyn.load("/Users/so11/phd/so11_nfs/somatic_ibd_p1/phylogenics/treemut.so")
  }
}
}

load_lib()


reconstruct_genotype_summary=function(phylo){
  dat=phylo$edge
  samples=phylo$tip.label
  N=length(samples)
  zeros=rep(0,N)
  profile=sapply(1:length(samples),function(i){tmp=zeros;tmp[i]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profile,edge_length=phylo$edge.length[1:N])
  #Create empty list where each element will correspond to a node
  muts=lapply(1:dim(dat)[1],function(i) c())
  map_node=match(1:max(dat[,2]),dat[,2])
  #Here we loop through the tips (samples) and get all the ancestral nodes for each tip.
  #Then for each node in the list we add the tip to that nodes list of tips, ultimately giving us a list of samples that share the node. 
  for(i in 1:N){
    parents=get_ancestral_nodes(i,edge = dat,exclude_root=TRUE)
    for(j in parents){
      muts[[map_node[j]]]=append(muts[[map_node[j]]],i)
    }
  }
  profiles=sapply(muts,function(x){tmp=zeros;tmp[x]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profiles,edge_length=phylo$edge.length,stringsAsFactors=FALSE)
  df=add_derived_profile_info(df,phylo$tip.label)
  list(df=df,samples=phylo$tip.label)
}

get_ancestral_nodes=function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}

add_derived_profile_info=function(profile_df,samples=sprintf("s%s",0:(nchar(profile_df$profile[1])-1))){
  base=rep(0,nchar(profile_df$profile[1]))
  samples_private=sapply(1:length(base),function(i){this=base;this[i]=1;paste(this,collapse="")})
  missing_samples=setdiff(samples_private,profile_df$profile)
  if(length(missing_samples)>0){
    profile_df=rbind(profile_df,data.frame(profile=missing_samples,edge_length=0))
  }
  
  profile_df$mut_count=nchar(profile_df$profile)-nchar(gsub("1","",profile_df$profile))
  profile_df$profile_int=lapply(profile_df$profile,split_profile)
  profile_df$id=1:dim(profile_df)[1]
  profile_df$label=profile_df$profile
  idx=which(profile_df$mut_count==1)
  profile_df$label[idx]=sapply(idx,function(i) samples[which(profile_df$profile_int[[i]]==1)])
  profile_df
}


assign_to_tree=function(mtr,dep,df,error_rate=rep(0.01,dim(mtr)[2]),maxits=5){
  if(!is.matrix(mtr)){
    nm=names(mtr)
    mtr=matrix(mtr,nrow=1)
    colnames(mtr)=nm
    dep=matrix(dep,nrow=1)
    colnames(dep)=nm
  }
  p.err=error_rate#c(error_rate,1e-10)
  p.err=p.err[match(df$samples,colnames(mtr))]
  mtr=mtr[,df$samples]
  dep=dep[,df$samples]
  if(!is.matrix(mtr)){
    nm=names(mtr)
    mtr=matrix(mtr,nrow=1)
    colnames(mtr)=nm
    dep=matrix(dep,nrow=1)
    colnames(dep)=nm
  }
  tree_genotypes=do.call("rbind",df$df$profile_int)
  el=rep(1,length(df$df$edge_length))
  loglik=rep(NA,maxits)
  n=dim(mtr)[1]
  for(i in 1:maxits){
    ol=el
    lik=get_likelihood_mtr_C(mtr,dep,tree_genotypes,el,p.error = p.err)
    edge_ml=apply(lik,1,which.max)
    #n=length(edge_ml)
    p=exp(lik-lik[(edge_ml-1)*n+1:n])##Substract max to just control range of lik... comes out in the wash later.
    p=p/rowSums(p)
    loglik[i]=sum(p*lik)
    el=colSums(p,na.rm=T)
    epsilon=sum(abs(el-ol))/dim(mtr)[1]
    cat("delta edge length=",epsilon,"\n")
    cat("Loglik=",loglik[i],"\n")
    if(epsilon<0.01/dim(mtr)[1]){
      break
    }
  }
  
  df$df$expected_edge_length=el
  df$df$edge_length=sapply(1:length(df$df$edge_length),function(i){length(which(edge_ml==i))})
  
  p_else_where=1-p[(edge_ml-1)*n+1:n]
  cat("Finished assigning mutations\ncalculating pvalues\n")
  pval=rep(NA,length(edge_ml))
  for(i in 1:length(pval)){
    pval[i]= get_mutation_assignment_pval(df$df,edge_ml[i],mtr[i,],dep[i,],p.err)
    if(i %% 1000 == 0){
      cat("On",i," of ",length(pval),"\n")
    }
  }
  
  list(df=df,lik=lik,summary=data.frame(edge_ml=edge_ml,pval=pval,p_else_where=p_else_where),p=p)
}



get_likelihood_mtr_C=function(mtr,depth,geno,el,p.error=rep(0.01,dim(mtr)[2])){
  if(dim(mtr)[2]!=dim(geno)[2]){
    stop("error: dimension mismatch betwee tree genotypes and mtr")
  }
  nmuts=dim(mtr)[1]
  nsamp=dim(mtr)[2]
  nbranch=dim(geno)[1]
  res = .C("likelihood_mtr", 
           mtr=as.integer(mtr),
           depth=as.integer(depth),
           geno=as.integer(geno),
           el=as.double(el),
           p.error=as.double(p.error),
           nmuts=as.integer(nmuts),
           nsamp=as.integer(nsamp),
           nbranch=as.integer(nbranch),
           lik=double(nmuts*nbranch)
  )
  matrix(res[["lik"]],ncol=nbranch)
}

get_mutation_assignment_pval=function(df,i,mtr,dep,p.error){
  profile_int=df$profile_int[[i]]
  idx.inside=which(profile_int==1)
  idx.outside=which(profile_int==0)
  V=0.5-p.error
  
  ##Calculate probability for observing <= total MTR at variant sites give VAF=0.5.  This is complicated by allowing different error rates for each sample (otherwise just pbinom)
  if(length(idx.inside)>0){
    ##We breakdown errors into unique categories to make calculation more tractable - could further discretise.
    v=unique(V[idx.inside])
    if(length(v)>4){
      stop("Too many error categories in get_pval: not tractable!")
    }
    mtri=mtr[idx.inside]
    depi=dep[idx.inside]
    vi=V[idx.inside]
    mtrii=sapply(v,function(v) sum(mtri[vi==v]))
    depii=sapply(v,function(v) sum(depi[vi==v]))
    p1=cpv(sum(mtrii),depii,v,lower.tail = TRUE)
    d1=1
  }else{
    p1=1
    d1=0
  }
  ##Calculate probability for observing >= total MTR at non-variant sites give VAF=0.  This is complicated by allowing different error rates for each sample (otherwise just pbinom)
  if(length(idx.outside)>0){
    v=unique(p.error[idx.outside])
    if(length(v)>4){
      stop("Too many error categories in get_pval: not tractable!")
    }
    mtri=mtr[idx.outside]
    depi=dep[idx.outside]
    vi=p.error[idx.outside]
    mtrii=sapply(v,function(v) sum(mtri[vi==v]))
    depii=sapply(v,function(v) sum(depi[vi==v]))
    p2=cpv(sum(mtrii),depii,v,lower.tail = FALSE)
    d2=1
  }else{
    p2=1
    d2=0
  }
  ##We are interested in whether either of the tests fails so we combine using conservative bonferroni.  
  ##We are performing 2 tests.  same as p.adjust(c(p1,p2),method = "bonferroni")
  pv=min((d1+d2)*min(c(p1,p2)),1)
  pv
}


cpv=function(mtrtot,depth,probs,lower.tail){
  ##Calculates the probability that there are a total of MTR mutant reads across N bins each with depth d_i and mutant prob p_i
  ##Does this by explicitly calculating  P(M | {d_i},{p_i)=sum_m P(mtr_i=m | d_i, p_i)*P(M-m | )
  n=length(depth)
  idx=order(depth)
  probs=probs[idx]
  depth=depth[idx]
  flag=ifelse(lower.tail,1,0)
  resk=list(mtrtot,probs,depth,n,lower.tail)
  res=.C("cumulate_binomial_pval",as.integer(round(mtrtot)),as.integer(round(depth)),as.double(probs),as.integer(n),as.integer(lower.tail),p=double(1))
  res[["p"]]
}

split_profile=function(profile){
  as.integer(unlist(strsplit(profile,split="")))
}


#####Crude Simulation code
simulate_reads_from_tree=function(nsamples,avgdepth,n_artifacts=0,p.error=0.01){
  if(length(p.error)>1){
    stop("p.error must be scalar here")
  }
  ##The following creates a random tree  and adds zero mutation outgroup (zeros)
  tree=rtree(nsamples)
  tree=bind.tree(tree,read.tree(text="(zeros:0);"))
  tree$edge.length=ceiling(200*tree$edge.length) ##Scales the #Muts per edge to the scale 0-200
  
  df=reconstruct_genotype_summary(tree)
  N=length(df$df$edge_length)
  dat=lapply(1:N,function(i) get_simulated_reads(df$df$profile_int[[i]],df$df$edge_length[i],avgdepth,p.error=p.error))
  ml=do.call("c",lapply(1:N,function(i) rep(i,df$df$edge_length[i])))
  if(n_artifacts>0){
    ##Add in some artefacts that aren't simulated according to the tree topology.
    NART=n_artifacts
    dat2=lapply(1:NART,function(i) get_simulated_reads(ifelse(runif(nsamples+1)<0.3,1,0),1,avgdepth,p.error=p.error))
    dat=c(dat,dat2)
    ml=c(ml,rep(NA,NART))
  }
  mtr=do.call("rbind",lapply(dat,function(x) x$mtr))
  depth=do.call("rbind",lapply(dat,function(x) x$depth))
  colnames(mtr)=df$samples
  colnames(depth)=df$samples
  p.error=rep(p.error,length(df$samples))
  idx.zeros=match("zeros",df$samples)
  mtr[,idx.zeros]=0
  depth[,idx.zeros]=10
  p.error[idx.zeros]=1e-6 ##This will force the inference code to keep the zeros branch as a zero length outgroup. 
  list(tree=tree,mtr=mtr,depth=depth,edge=ml,df=df,p.error=p.error)
}

get_simulated_reads=function(geno,n,depth,p.error=0.01){
  if(n==0){
    emptymat=matrix(1,nrow=0,ncol=length(geno))
    return(list(mtr=emptymat,depth=emptymat))
  }
  geno=rep(geno,n)
  depth=rpois(length(geno),depth)
  mtr=rbinom(length(depth),depth,prob=ifelse(geno==0,p.error,0.5-p.error))
  list(mtr=matrix(mtr,nrow=n,byrow = TRUE),
       depth=matrix(depth,nrow=n,byrow = TRUE)
  )
}

run_sim=function(nsamp,depth){
  dat=simulate_reads_from_tree(nsamp,depth)
  res=assign_to_tree(dat$mtr,dat$depth,dat$df,error_rate=dat$p.error)##
  list(edge_length_orig=dat$df$df$edge_length,
       edge_length_inferred=res$df$df$edge_length,
       expected_edge_length_inferred=res$df$df$expected_edge_length,
       edge_idx_orig=dat$edge,
       edge_idx_ml=res$summary$edge_ml)
}

plot_sim_result=function(sim){

  
  par(mfrow=c(2,2))
  mismatch_prop=length(which(sim$edge_idx_orig!=sim$edge_idx_ml))/length(sim$edge_idx_orig)
  plot(sim$edge_length_orig,sim$edge_length_inferred,xlab="Original Edge Length",ylab="Inferred Edge Length",
       main=sprintf("Hard Assigned edge length vs Orig (sd=%3.2f)\n Mismatch Proportion=%5.4f",
                    sd(sim$edge_length_orig-sim$edge_length_inferred),mismatch_prop)
  )
  abline(a=0,b=1)
  plot(sim$edge_length_orig,sim$expected_edge_length_inferred,xlab="Original Edge Length",ylab="Expected Edge Length",
       main=sprintf("Expected edge length (Soft Assigned) \nvs Orig (sd=%3.2f)",sd(sim$edge_length_orig-sim$expected_edge_length_inferred))
  )
  abline(a=0,b=1)
  simdf=data.frame(deviation=sim$expected_edge_length_inferred-sim$edge_length_orig,
                   edge_length_orig=sim$edge_length_orig)
  simdf=simdf[order(simdf$edge_length_orig),]
  
  with(simdf,plot(edge_length_orig,deviation,
       xlab="Original Edge Length",ylab="Deviation",
       main="Deviation vs Original Edge Length"))
  ##Add LOESS to see if there is any systematic patterns
  loess_fit <- with(simdf,loess(deviation ~ edge_length_orig))
  pl=predict(loess_fit,se = T)
 
  plu=pl$fit+1.96*pl$se
  pll=pl$fit-1.96*pl$se
  
  
  shade_between(simdf$edge_length_orig,plu,pll,adjustcolor("blue",0.3))
  with(simdf,lines(edge_length_orig, pl$fit, col = "blue",lwd=2))
  
  plot(simdf$edge_length_orig,abs(simdf$deviation),
       xlab="Original Edge Length",ylab="Absolute Deviation",
       main="Absolute Deviation vs Original Edge Length")
  
  
  
}

shade_between=function(x,y1,y2,color){
  polygon(c(x, rev(x), x[1]), c(y1, rev(y2), y1[1]), 
          col = color,border = color) 
}



#-------------------------------------------------
# Functions for sensitivity adjustment of the branhces
#-------------------------------------------------


## Input: tree is a phylo object with mutation counts listed in the branch.length attribute.
##        coverage is a two-column dataframe where the first column contains sample names that match the tip labels
##        of the tree and the second column is the coverage of that sample.
##        vaf is either a two-column dataframe where the first column contains sample names that match the tip labels
##        of the tree and the second column is the median vaf of that sample, a numeric vector
##        containing vafs for the samples in coverage (in the same order) or a single number.
## Output: A phylo object with adjusted mutation counts listed in the branch.length attribute.

## This has really been depricated now. Using the truncated
## binomial correction as suggested by Inigo. 
adjust_branch_lengths <- function(tree, coverage, vaf=0.5, method=c("betaBinom", "Poisson", "ml"), read_threshold= 4) {
  tree_df <- fortify(tree)
  
  colnames(coverage) <- c("label", "Coverage")
  if(is.data.frame(vaf)) {
    colnames(vaf) <- c("label", "vaf")
    coverage <- merge(coverage, vaf)
  } else {
    coverage$vaf <- vaf
  }
  # sometimes the matched normal has no mutations. Just dealing with that. 
  coverage$vaf[is.na(coverage$vaf)] <- 0.5
  ## match the coverage with the tip labels.
  y <- merge(tree_df, coverage, all.x=T)
  y <- y[order(y$node),]
  
  
  ## For each node that isn't a leaf, set the coverage for that branch as the sum of the coverage of
  ## all it's children. Set the vaf for the branch as the average median vaf for the children.
  if(nrow(tree_df) > 3) {
    for(i in (sum(tree_df$isTip)+2):nrow(tree_df)) {
      y$Coverage[i] <- sum(y$Coverage[find_children(i, tree_df)])
      y$vaf[i] <- mean(y$vaf[find_children(i, tree=tree_df)])
    }
  } 
  
  # Set the coverage and vaf of the root (not really needed since it will always have 0 mutations)
  y$Coverage[(sum(tree_df$isTip)+1)] <- sum(coverage$Coverage)
  y$vaf[(sum(tree_df$isTip)+1)] <- 0.5
  
  ## Use the selected method to estimate the sensitivity
  sensitivity <- numeric()
  for(i in 1:nrow(tree_df)) {
    
    if(method=="betaBinom") {
      ## When rho=0, this reduces to binomial distribution. Have the coverage follow a Poisson distribution with lambda as mean cov. 
      sensitivity[i] <- sum(rbetabinom(10000, size=round(y$Coverage[i]), prob=y$vaf[i], rho=0) >=read_threshold)/10000
    } else if(method=="Poisson") {
      ## Use half the coverage since mutations are heterozygous
      ## How often do we see at least read_threshold reads supporting the variant?
      sensitivity[i] <- sum(rpois(10000, lambda=y$Coverage[i]*y$vaf[i]) >=read_threshold)/10000
    } else if(method=="ml") {
      stop("Maximum likelihood sensitivity correction not implemented yet")
      ## will need a linear model that predicts sensitivity from coverage and vaf.
      ## Not sure Simon's approach is applicable for my scenario as crypts don't represent
      ## two samples from the same clone. Given a pair of crypts you expect both to have their
      ## own set of unique mutations. Simon's approach assumes you've sampled exactly the same
      ## clone twice. You punish for unique mutations.
    } else {
      stop(paste("The method", method, "is not implemented for this function."))
    }
  }
  
  ## Adjust the branch.length attribute for sensitivity
  tree$edge.length <- tree_df$branch.length[tree$edge[,2]]/sensitivity[tree$edge[,2]]
  
  ## Return updated tree
  return(tree)
}

find_children = function(node, tree=tree_df){
  # Function to find terminal tips belonging to 
  # branch in phylogenetic tree 
  child_nodes = tree$node[tree$parent==node]
  tips = c()
  for (k in 1:length(child_nodes)){
    if (tree$isTip[tree$node==child_nodes[k]]){
      tips=c(tips,child_nodes[k])
    }
    else{
      tips=c(tips,find_children(child_nodes[k])) 
    }
  }
  return(tips)
}



### Functions for sensitivity correction by fitting a
### truncated beta-binomial. 

trunc_binom_adjustment <- function(tree, NR, NV, coverage, hard_limit=3) {
  tree_df <- fortify(tree)
  
  sensitivities <- numeric()
  for(i in 1:ncol(NR)) {
    
    ## If there's a matched normal set the sensitivity for it
    ## as 1. 
    if(grepl("No", colnames(NR))[i]) {
      sensitivities <- c(sensitivities, 1)
      next
    }
    
    ## Only use variants with 4 or more reads
    NR_4=NR[NV[,i]>hard_limit,i]
    NV_4=NV[NV[,i]>hard_limit,i]
    
    res = binom_mix(NV_4,NR_4,mode="Truncated") 
    
    sensitivities <- c(sensitivities, mean(unlist(lapply(rpois(n=100000,lambda=cov[i,"Coverage"]),function(x) pbinom(q=hard_limit,size=x,p=res$p,lower.tail = F)))))
  }
  
  sens <- data.frame(label=colnames(NR), sensitivity=sensitivities)
  
  ## match the sensitivity with the tip labels.
  y <- merge(tree_df, sens, all.x=T)
  y <- y[order(y$node),]
  
  ## Set the sensitivity of the root as 1.
  y$sensitivity[(sum(tree_df$isTip)+1)] <- 1
  
  ## For each node that isn't a leaf, set the coverage for that branch as the sum of the coverage of
  ## all it's children. Set the vaf for the branch as the average median vaf for the children.
  if(nrow(tree_df) > 3) {
    for(i in (sum(tree_df$isTip)+2):nrow(tree_df)) {
      child_sens <- y$sensitivity[find_children(i, tree_df)]
      inv_child_sens <- 1-child_sens
      y$sensitivity[i] <- 1-prod(inv_child_sens)
    }
  } 
  
  ## Adjust the branch.length attribute for sensitivity
  tree$edge.length <- tree_df$branch.length[tree$edge[,2]]/y$sensitivity[tree$edge[,2]]
  
  ## Return updated tree
  return(tree)
}

## Functions from Tim Coorens
##################################

## Define the truncated binomial distribution
dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
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
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
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

binom_mix = function(x,size,nrange=1:5,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
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
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
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
#-------------------------------------------------









