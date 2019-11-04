
## This script does the following:
## 1. Read in a MPBoot tree, genotype matrix and the nr of reads supporting mutation/ref alleles.
## 2. Use maximum Likelihood to assing each mutation to a branch of the tree
## 3. Plot the tree after this
## 4. Adjust the branch lengths by VAF and coverage. Plot a new tree and sum the mutation counts for
## each sample to be used in lmm of mutation burden.
## 5. Assign mutations to branches and write this out for signature extraction. 



## Define parameters
if(file.exists("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt")) {
  lib_location="/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm3/"
  .libPaths( c(lib_location))
  library(ape)
  library(ggtree)
  library(VGAM)
#  .libPaths( c(lib_location, .libPaths()) )
  source("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/phylogenics/treemut.R")
  #biopsy_meta <- read.table("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/sample_info/biopsy_metadata.txt", h=T)
  crypt_meta <- read.table("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt", h=T)
  
  parameters <- commandArgs(TRUE)
  PID <- parameters[1]
  gt_matrix <- parameters[2]
  treeFile <- parameters[3]
  NR_file <- parameters[4]
  NV_file <- parameters[5]
  vaf_file <- parameters[6]
  output_dir <- parameters[7]
  mutType <- parameters[8]
  
} else {
  #biopsy_meta <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/sample_info/biopsy_metadata.txt", h=T)
  crypt_meta <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt", h=T)
  ## For testing: 
  PID <- "patient28"
  gt_matrix <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/qc/binomial_filters/patient28/patient28_passing_snps_genotype_binary_matrix.txt"
  treeFile <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/phylogenics/MPBoot/snv/consensus_trees/patient28/patient28.treefile"
  output_dir <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/phylogenics/MPBoot/snv/"
  NR_file <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/qc/binomial_filters/patient28/patient28_NR_pass.txt"
  NV_file <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/qc/binomial_filters/patient28/patient28_NV_pass.txt"
  vaf_file <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/qc/binomial_filters/patient28/patient28_medVAF_and_Unadj_MutCount.txt"
  source("/Users/so11/phd/so11_nfs/somatic_ibd_p1/phylogenics/treemut.R")
  library(ape)
  library(ggtree)
  library(VGAM)
}

options(stringsAsFactors = F)

## Read in files: 
NR <- read.table(NR_file, h=T)
NV <- read.table(NV_file, h=T)
Genotype <- read.table(gt_matrix, h=T)
tree=read.tree(treeFile)


## Clean up the matrices a bit
NR <- NR[complete.cases(NR), ]
rownames(NR) <- NR$mutID
NR$mutID <- NULL
colnames(NR) <- gsub("_DEP","",colnames(NR))

NV <- NV[complete.cases(NV), ]
rownames(NV) <- NV$mutID
NV$mutID <- NULL
colnames(NV) <- gsub("_DEP","",colnames(NV))

## Clean the tree as well. 
tree$tip.label=gsub("_VAF","",tree$tip.label)
tree$edge.length=rep(1,nrow(tree$edge))
tree=drop.tip(tree,"Ancestral")
NR_flt = as.matrix(NR[rownames(Genotype),tree$tip.label])
NV_flt = as.matrix(NV[rownames(Genotype),tree$tip.label])

## Get the meta-data sorted:
meta <- crypt_meta[, c("crypt_ID", "biopsy_ID", "location", "type")]



## Assign mutations to the tree
df=reconstruct_genotype_summary(tree)
res=assign_to_tree(df=df,
                   mtr=NV_flt,
                   dep=NR_flt)


## Make the edge lengths equal 
edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<0.01])
#edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<=1])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero

tree2=tree
tree2$edge.length=as.numeric(edge_length)

system(paste("mkdir -p ", output_dir, "consensus_trees/", PID, "/", sep=""))
system(paste("mkdir -p ", output_dir, "tree_plots/", PID, "/", sep=""))

write.tree(tree2, paste(output_dir, "consensus_trees/", PID, "/", PID, "_MutAssign_unadj.tree", sep=""))

png(paste(output_dir, "tree_plots/", PID, "/", PID, "_MutAssign_unadj.png", sep=""), width = 1000, height=600)
  plot <- ggtree(tree2) 
  plot <- plot %<+% meta + geom_tiplab(aes(color=paste(location, type)),vjust=-0.3, size=6) +
    geom_tippoint(aes(color=paste(location, type)), alpha=0.25)
  plot <- plot +theme_tree2()+xlim(0,max(fortify(tree2)$x)*1.3) + 
    theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size=16))
  print(plot)
dev.off()

png(paste(output_dir, "tree_plots/", PID, "/", PID, "_MutAssign_unadj_nodeNrs.png", sep=""), width = 1000, height=600)
plot <- ggtree(tree2) 
plot <- plot %<+% meta + geom_tiplab(aes(color=paste(location, type)),vjust=-0.3, size=6) +
  geom_tippoint(aes(color=paste(location, type)), alpha=0.25)
plot <- plot +theme_tree2()+xlim(0,max(fortify(tree2)$x)*1.3) + 
  geom_text2(aes(label=node), hjust=-.3, vjust=0.9, color="red") + 
  theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size=16))
print(plot)
dev.off()


## Next I'm going to adjust the mutation count for coverage and VAF of the crypt
## See functions in treemut.R

## Need the coverage and the median vaf:
cov <- crypt_meta[crypt_meta$crypt_ID %in% tree2$tip.label, c("crypt_ID", "Coverage")]
vafs <- read.table(vaf_file, h=T)
if(mutType=="snv") {
  vafs <- vafs[, c("crypt_ID", "snvMedianVAF")]
} else if(mutType=="indel") {
  vafs <- vafs[, c("crypt_ID", "indelMedianVAF")]
}

vafs$crypt_ID <- gsub("_VAF", "", vafs$crypt_ID)
vafs$mutation_counts <- NULL

tree_df <- fortify(tree2)
adj_tree <- adjust_branch_lengths(tree2, coverage = cov, vaf = vafs, method="Poisson")

write.tree(adj_tree, paste(output_dir, "consensus_trees/", PID, "/", PID, "_MutAssign_Poisson_adj.tree", sep=""))


png(paste(output_dir, "tree_plots/", PID, "/", PID, "_MutAssign_Poisson_adjusted.png", sep=""), width = 1000, height=600)
plot <- ggtree(adj_tree) 
plot <- plot %<+% meta + geom_tiplab(aes(color=paste(location, type)),vjust=-0.3, size=6) +
  geom_tippoint(aes(color=paste(location, type)), alpha=0.25)
plot <- plot +theme_tree2()+xlim(0,max(fortify(adj_tree)$x)*1.3) + 
  theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size=16))
print(plot)
dev.off()

adj_df <- fortify(adj_tree)

## For each node that is a leaf, sum up all branch lengths
adj_mutCount <- adj_df[adj_df$isTip, c("label","x")]
colnames(adj_mutCount) <- c("crypt_ID", "adj_mutCount")

#system(paste("mkdir -p ", output_dir, "adj_mutCounts/", PID, "/", sep=""))
write.table(adj_mutCount, paste(output_dir,"adj_mutCounts/", PID,"_Poisson_adjusted.txt", sep="" ), quote = F, row.names = F)

## Adjusting branch lengths using truncated binomials:
trunc_adj_tree <- trunc_binom_adjustment(tree=tree2,NR=NR, NV=NV, coverage=cov)

write.tree(trunc_adj_tree, paste(output_dir, "consensus_trees/", PID, "/", PID, "_MutAssign_truncBinom_adj.tree", sep=""))

png(paste(output_dir, "tree_plots/", PID, "/", PID, "_MutAssign_truncBinom_adjusted.png", sep=""), width = 1000, height=600)
plot <- ggtree(trunc_adj_tree) 
plot <- plot %<+% meta + geom_tiplab(aes(color=paste(location, type)),vjust=-0.3, size=6) +
  geom_tippoint(aes(color=paste(location, type)), alpha=0.25)
plot <- plot +theme_tree2()+xlim(0,max(fortify(trunc_adj_tree)$x)*1.3) + 
  theme(legend.position="right", legend.title = element_blank(), legend.text = element_text(size=16))
print(plot)
dev.off()

trunc_adj_df <- fortify(trunc_adj_tree)

## For each node that is a leaf, sum up all branch lengths
trunc_adj_mutCount <- trunc_adj_df[trunc_adj_df$isTip, c("label","x")]
colnames(trunc_adj_mutCount) <- c("crypt_ID", "trunc_adj_mutCount")

write.table(trunc_adj_mutCount, paste(output_dir,"adj_mutCounts/", PID,"_truncBinom_adjusted.txt", sep="" ), quote = F, row.names = F)



## Write out a file assigning mutations to branches. 
system(paste("mkdir -p ", output_dir, "branch_mut_assignment/", PID, "/", sep=""))

con = textConnection(rownames(Genotype))
Mutations_per_branch = read.table(con,sep="_")
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch$Branch[rowSums(Genotype>0.3)==ncol(Genotype)]=0 #Root
Mutations_per_branch=Mutations_per_branch[Mutations_per_branch$Branch!=0&res$summary$p_else_where<0.01,]
Mutations_per_branch$Patient = PID
Mutations_per_branch$SampleID = paste(PID,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste(output_dir, "branch_mut_assignment/", PID, "/",PID, "_ML_branch_assignment.txt", sep=""),quote=F,row.names=F,sep="\t")



