
## The purpose of this script is to perform HDP-signature extraction for substitutions

## See Henry's directory here: /lustre/scratch117/casm/team154/hl11/colon_temp/signatures_May_2018/hdp/substitutions_2018.04.30/

if(file.exists("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt")) {
  lib_location="/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm3/"
  .libPaths( c(lib_location, .libPaths()) )
  library(hdp, lib.loc=lib_location)
  output_dir <- "/lustre/scratch114/projects/crohns/somatic_ibd_p1/signature_extraction/hdp/snvs/"
  # read in PCAWG sigs
  all_sigs <- read.csv("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/resources/pcawk_signatures/sigProfiler_SBS_signatures_summing_to_one.csv", header = T, row.names = 1, stringsAsFactors = F)
} else {
  output_dir <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/"
  # read in PCAWG sigs
  all_sigs <- read.csv("/Users/so11/phd/so11_nfs/somatic_ibd_p1/resources/pcawk_signatures/sigProfiler_SBS_signatures_summing_to_one.csv", header = T, row.names = 1, stringsAsFactors = F)
  
  library(hdp)
  # see package vignettes with
  # browseVignettes('hdp')
}


load(paste(output_dir, "branch_mutational_matrix.txt", sep=""))

jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
randseed <- 42*(jobIndex) + 1
print(paste0("Job index:", jobIndex))
print(paste0("Random seed:", randseed))


tnt <- mut_mat_t
newColNames <- paste(substr(colnames(mut_mat_t), 3,3), ".", substr(colnames(mut_mat_t), 5,5), ".in.", substr(colnames(mut_mat_t), 1,1), substr(colnames(mut_mat_t), 3,3), substr(colnames(mut_mat_t), 7,7), sep="")
colnames(tnt) <- newColNames
patients <- sapply(strsplit(rownames(tnt), "_"), "[[", 1)


## Add Henry's data to the mix

if(FALSE) {
  #henry_sbs <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/henry_files/sbs_category_counts.txt", h=T) 
  henry_sbs <- read.table("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/henry_files/sbs_category_counts.txt", h=T)
  tnt <- rbind(tnt, henry_sbs)
  henry_patients <- sapply(strsplit(rownames(henry_sbs), ":"), "[[", 1)
  patients <- c(patients, henry_patients)
  
}



## Filter the signatures.
# I only want to keep a subset of them. I am picking the ones that have been found in colorectal cancer, as documented in the PCAWG SBS vignette 2018.04.27, downloaded from synapse 2018.04.30
#gdsigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS6", "SBS9", "SBS10a", "SBS10b", "SBS13", "SBS15", "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS20", "SBS21", "SBS25", "SBS26", "SBS28", "SBS30", "SBS37", "SBS40", "SBS41", "SBS43", "SBS44", "SBS45", "SBS46", "SBS49") 
gdsigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS25", "SBS28","SBS32", "SBS30", "SBS37", "SBS40", "SBS41", "SBS43", "SBS45", "SBS49")
sigs <- as.matrix(all_sigs[,gdsigs])

prior_pseudoc <- rep(1000, ncol(sigs)) # this is meant to be a vector of pseudocounts contributed by each prior distribution.
test_prior <- hdp_prior_init(sigs, prior_pseudoc, hh=rep(1, 96), alphaa=c(1, 1), alphab=c(1, 1))

# add three more CPs for the data that we will add: 
# one for the branch from the grandparent to the MRCA of the patients, 
# one for all the branches from the MRCA of the patients to each patient
# one to go from each patient to all their child nodes (same cp for all patients)
test_prior <- hdp_addconparam(test_prior, alphaa=c(1,1,1), alphab=c(1,1,1))



# now add the hierarchical structure of the dataset.
# I want: 
# ppindex 0 for the top node. - ALREADY PRESENT
# ppindex 1 for each of the 28 signatures. - ALREADY PRESENT
# ppindex 1 for the MRCA of all the patients. - TO ADD
# ppindex 30 (i.e. the MRCA of all the patients) for all 42 patients. - TO ADD
# and then each branch has the ppindex of the patient that it comes from. - TO ADD

ptsmrcappindex <- 1
ptsppindices <- rep(c(1 + ncol(sigs) + 1), length(unique(patients)))

branchppindices <- as.numeric(unlist(sapply(unique(patients), function(patient) {
  tnum <- length(which(patients==patient))
  tindex <- (1 + ncol(sigs) + 1) + which(unique(patients)==patient)
  rep(tindex, tnum)
})))

newppindices <- c(ptsmrcappindex, ptsppindices, branchppindices)

# For concentration parameters, give one for every level of the hierarchy:
# cpindex 1 for the grandparent. - ALREADY PRESENT
# cpindex 2 for all the signatures - ALREADY PRESENT
# cpindex 3 for the parent of the patients - NEED TO ADD
# cpindex 4 for all the patients - NEED TO ADD
# cpindex 5 for all the patients to their branches - NEED TO ADD
newcpindices <- c(3, rep(4, length(unique(patients))), rep(5, nrow(tnt)))

# add dp nodes: 
# one as the parent of all the patients
# one for each of the 42 patients
# one for every single branch.
# i.e. this should be the same as the number of new ppindicies and cpindices
test_prior <- hdp_adddp(test_prior,
                        numdp = length(newppindices),
                        ppindex = newppindices,
                        cpindex = newcpindices)

# need to make sure that I am adding data in a way that matches the data to the terminal nodes.
# dpindices are the indices of the terminal nodes. The first terminal node is number 73.
test_prior <- hdp_setdata(hdp=test_prior, dpindex=(max(newppindices) + 1:nrow(tnt)), data=tnt)


# run chain
test_pr <- dp_activate(test_prior, dpindex=(1+ncol(sigs)+1):numdp(test_prior), initcc=(ncol(sigs)+round(ncol(sigs)/10)))
test_chlist <- hdp_posterior(test_pr, burnin=100000, n=100, space=2000, cpiter=3, seed=randseed)

# save data.
assign(paste0("hdp_", jobIndex), test_chlist)

save(list=paste0("hdp_", jobIndex), file=paste(output_dir, "hdpout_", jobIndex, ".RData", sep=""))








