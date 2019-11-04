
## The purpose of this script is to run binomial filters to remove germline variants from somatic mutation calls.


  parameters <- commandArgs(TRUE)
  ## Crypt fission rate in the unit fissions per crypt per year.
  snp_file <- parameters[1]
  indel_file <- parameters[2]
  annotation_file <- parameters[3]
  sex <- parameters[4]
  cut_off <- as.numeric(parameters[5])
  nr_crypts_threshold <- as.numeric(parameters[6])
  script_dir <- parameters[7]
  output_dir <- parameters[8]
  patient <- parameters[9]


if(FALSE) {
patient="patient45"
snp_file=paste("/Users/so11/phd/so11_lustre/somatic_ibd_p1/pileups/ibd_cohort/", patient, "/snp/", patient, "_snp_binomial_input.txt", sep="")
indel_file=paste("/Users/so11/phd/so11_lustre/somatic_ibd_p1/pileups/ibd_cohort/", patient, "/indel/", patient, "_indel_binomial_input.txt", sep="")
annotation_file=paste("/Users/so11/phd/so11_lustre/somatic_ibd_p1/pileups/ibd_cohort/", patient, "/", patient, "_allinfo_combined.txt", sep="")

sex="female"
cut_off <- 3
nr_crypts_threshold <- 5
script_dir="/Users/so11/phd/so11_nfs/somatic_ibd_p1/qc/binomial_filters/"
output_dir=paste("/Users/so11/phd/so11_lustre/somatic_ibd_p1/qc/binomial_filters/", patient, "/", sep="") 

}
  
## When there are too few crypts the cut-off for the exact binomial needs to be low so I don't go ahead and filter everything.

  
source(paste(script_dir, "binomial_filters_function_archive.R", sep=""))


## Hard-coded parameters
gt_matrix_lower <- 0.05
gt_matrix_upper <- 0.2
options(stringsAsFactors = F)

snps <- read.table(snp_file, h=T,stringsAsFactors=F)
indels <- read.table(indel_file, h=T, stringsAsFactors=F)

nrCrypts <- (ncol(snps) - 15)/15
if(nrCrypts<nr_crypts_threshold) {
  cut_off <- 2
}

depth_threshold <- 5
mutant_reads_threshold <- 3

snps <- snps[, !grepl("P42B1_3", colnames(snps))]
indels <- indels[, !grepl("P42B1_3", colnames(indels))]
snps <- snps[, !grepl("P42B1_4", colnames(snps))]
indels <- indels[, !grepl("P42B1_4", colnames(indels))]

snps <- snps[, !grepl("P30B2", colnames(snps))]
indels <- indels[, !grepl("P30B2", colnames(indels))]

snps <- snps[, !grepl("P54B2b_9", colnames(snps))]
indels <- indels[, !grepl("P54B2b_9", colnames(indels))]
snps <- snps[, !grepl("P54B2b_6", colnames(snps))]
indels <- indels[, !grepl("P54B2b_6", colnames(indels))]

snps <- snps[, !grepl("P10B2_23", colnames(snps))]
indels <- indels[, !grepl("P10B2_23", colnames(indels))]
snps <- snps[, !grepl("P10B2_22", colnames(snps))]
indels <- indels[, !grepl("P10B2_22", colnames(indels))]
snps <- snps[, !grepl("P10B2_24", colnames(snps))]
indels <- indels[, !grepl("P10B2_24", colnames(indels))]

snps <- snps[, !grepl("P22B3", colnames(snps))]
indels <- indels[, !grepl("P22B3", colnames(indels))]
#snps <- snps[, !grepl("P36B2a_4", colnames(snps))]
#indels <- indels[, !grepl("P36B2a_4", colnames(indels))]
snps <- snps[, !grepl("P48B1", colnames(snps))]
indels <- indels[, !grepl("P48B1", colnames(indels))]

snps <- snps[, !grepl("P48B2a_1", colnames(snps))]
indels <- indels[, !grepl("P48B2a_1", colnames(indels))]

snps <- snps[, !grepl("P46B2_1", colnames(snps))]
indels <- indels[, !grepl("P46B2_1", colnames(indels))]
snps <- snps[, !grepl("P46B2_3", colnames(snps))]
indels <- indels[, !grepl("P46B2_3", colnames(indels))]


NR_snp = snps[,grepl("DEP",colnames(snps))]
WTR_snp= snps[,grepl("WTR",colnames(snps))]

NR_indels = indels[,grepl("DEP",colnames(indels))]
WTR_indels= indels[,grepl("WTR",colnames(indels))]

## A bug was introduced in the pileup step as I added a "matched normal" for samples that weren't run with
## matched normals. I added it for BRASS filtering but it causes the same sample to be included twice in the 
## pileup. 
cryptNames_indel <- gsub("_DEP.[0-9]$", "", colnames(NR_indels))
cryptNames_indel <- gsub("_DEP", "", cryptNames_indel)
colnames(NR_indels) <- gsub("_DEP", "", colnames(NR_indels))
NR_indels <- NR_indels[, colnames(NR_indels) %in% cryptNames_indel]
colnames(WTR_indels) <- gsub("_WTR", "", colnames(WTR_indels))
WTR_indels <- WTR_indels[, colnames(WTR_indels) %in% cryptNames_indel]


cryptNames <- gsub("_DEP.[0-9]$", "", colnames(NR_snp))
cryptNames <- gsub("_DEP", "", cryptNames)
colnames(NR_snp) <- gsub("_DEP", "", colnames(NR_snp))
NR_snp <- NR_snp[, colnames(NR_snp) %in% cryptNames]
colnames(WTR_snp) <- gsub("_WTR", "", colnames(WTR_snp))
WTR_snp <- WTR_snp[, colnames(WTR_snp) %in% cryptNames]

NV_snp <- NR_snp-WTR_snp
rownames(NV_snp)=rownames(NR_snp)=rownames(WTR_snp)=c(paste(snps$Chrom, snps$Pos, snps$Ref, snps$Alt, sep="_"))

NV_indels <- NR_indels-WTR_indels
rownames(NV_indels)=rownames(NR_indels)=rownames(WTR_indels)=c(paste(indels$Chrom, indels$Pos, indels$Ref, indels$Alt, sep="_"))



NR <- rbind(NR_snp, NR_indels)
WTR <- rbind(WTR_snp, WTR_indels)
NV <- rbind(NV_snp, NV_indels)

vaf_indels <- indels[,grepl("VAF",colnames(indels))]
colnames(vaf_indels ) <- gsub("_VAF", "", colnames(vaf_indels ))
vaf_indels <- vaf_indels[, colnames(vaf_indels) %in% cryptNames_indel]

vaf_snp <- snps[,grepl("VAF",colnames(snps))]
colnames(vaf_snp) <- gsub("_VAF", "", colnames(vaf_snp))
vaf_snp <- vaf_snp[, colnames(vaf_snp) %in% cryptNames]

VAFs <- rbind(vaf_snp, vaf_indels)
rownames(VAFs) <- rownames(NR)

depth <- NR>=depth_threshold
mutant_reads <- NV>=mutant_reads_threshold
both <- depth & mutant_reads

## A function to filter out germline variants.
## Does very poorly when the aggregate coverage is low (few samples per patient)
germline_exact=exact.binomial(NV=NV,NR=NR,gender=sex,cutoff=-cut_off)
germline_exact[is.na(germline_exact)] <- FALSE
germline_binom <- beta.binom.filter(NR, NV)
germline_binom[is.na(germline_binom)] <- FALSE

germline <- !germline_exact & germline_binom

both[!germline, ] <- FALSE
x <- VAFs
for(i in 1:ncol(VAFs)) {
  for(j in 1:nrow(VAFs)) {
    if(both[j,i]==FALSE) {
      x[j,i] <- NA
    }
  }
}

passing_vafs <- VAFs[germline,]
failing_vafs <- VAFs[!germline,]
should_be_counted <- x[germline,]

NR_pass <- NR[germline,]
NV_pass <- NV[germline,]
NR_pass$mutID <- rownames(NR_pass)
NV_pass$mutID <- rownames(NV_pass)

write.table(NR_pass, paste(output_dir, patient, "_NR_pass.txt", sep=""), quote=F, row.names = F)
write.table(NV_pass, paste(output_dir, patient, "_NV_pass.txt", sep=""), quote=F, row.names = F)

passing_vafs <- passing_vafs[complete.cases(passing_vafs), ]
failing_vafs <- failing_vafs[complete.cases(failing_vafs), ]

global_vafs <- numeric()
for(j in 1:nrow(NR)) {
  global_vafs[j] <- sum(NV[j,])/sum(NR[j,])
}
names(global_vafs) <- rownames(NR)

global_vafs_pass <- global_vafs[germline]
global_vafs_fail <- global_vafs[!germline]
global_vafs_pass <- global_vafs_pass[!is.na(global_vafs_pass)]
global_vafs_fail <- global_vafs_fail[!is.na(global_vafs_fail)]

global_vafs_pass_snps <- subset(global_vafs_pass, names(global_vafs_pass) %in% rownames(NR_snp))
global_vafs_pass_indel <- subset(global_vafs_pass, names(global_vafs_pass) %in% rownames(NR_indels))
global_vafs_pass_snps <- global_vafs_pass_snps[!is.na(global_vafs_pass_snps)]
global_vafs_pass_indel <- global_vafs_pass_indel[!is.na(global_vafs_pass_indel)]

# This is really the average VAF, not the global vaf. 
#global_vafs_pass <- apply(passing_vafs, 1, sum) / ncol(passing_vafs)
#global_vafs_fail <- apply(failing_vafs, 1, sum) / ncol(failing_vafs)

g1 <- ggplot(data.frame(global_vafs_pass), aes(x=global_vafs_pass)) + geom_histogram()
g2 <- ggplot(data.frame(global_vafs_fail), aes(x=global_vafs_fail)) + geom_histogram()
g3 <- ggplot(data.frame(global_vafs_pass_snps), aes(x=global_vafs_pass_snps)) + geom_histogram()
g4 <- ggplot(data.frame(global_vafs_pass_indel), aes(x=global_vafs_pass_indel)) + geom_histogram()

png(paste(output_dir, patient, "_global_vaf_pass_all.png", sep=""))
  print(g1)
dev.off()

png(paste(output_dir, patient, "_global_vaf_fail_all.png", sep=""))
print(g2)
dev.off()

png(paste(output_dir, patient, "_global_vaf_pass_snps.png", sep=""))
print(g3)
dev.off()

png(paste(output_dir, patient, "_global_vaf_pass_indels.png", sep=""))
print(g4)
dev.off()

global_vafs_df <- data.frame(ID=c(rownames(NR)[germline], rownames(NR)[!germline]),
                                  global_VAF=c(global_vafs_pass, global_vafs_fail),
                                  status=c(rep("PASS", length(global_vafs_pass)), rep("FAIL", length(global_vafs_fail))))

write.table(global_vafs_df, paste(output_dir, patient, "_global_vafs.txt", sep=""), quote=F, row.names = F)

passing_vafs_snp <- subset(passing_vafs, rownames(passing_vafs) %in% rownames(NR_snp))
passing_vafs_indels <- subset(passing_vafs, rownames(passing_vafs) %in% rownames(NR_indels))

should_be_counted_snp <- subset(should_be_counted, rownames(should_be_counted) %in% rownames(NR_snp))
should_be_counted_indels <- subset(should_be_counted, rownames(should_be_counted) %in% rownames(NR_indels))

median_vafs_total <- apply(should_be_counted, 2, find_median_vaf, lowest=gt_matrix_upper)
mutation_counts_total <- apply(should_be_counted, 2, Mutcount, lower=gt_matrix_upper)

median_vafs_snp <- apply(should_be_counted_snp, 2, find_median_vaf, lowest=gt_matrix_upper)
mutation_counts_snp <- apply(should_be_counted_snp, 2, Mutcount, lower=gt_matrix_upper)

median_vafs_indels <- apply(should_be_counted_indels, 2, find_median_vaf, lowest=gt_matrix_upper)
mutation_counts_indels <- apply(should_be_counted_indels, 2, Mutcount, lower=gt_matrix_upper)

summary <- data.frame(crypt_ID=names(median_vafs_total), totalMedianVAF=median_vafs_total, totalMutCount=mutation_counts_total,
                      snvMedianVAF=median_vafs_snp, snvMutCount=mutation_counts_snp, indelMedianVAF=median_vafs_indels, indelMutCount=mutation_counts_indels)

write.table(summary, paste(output_dir, patient, "_medVAF_and_Unadj_MutCount.txt", sep=""), quote=F, row.names = F)



## Read in the Annotation files from Mathijs and write back out the subst of patients that pass.
## The annotation file only covers subs
if(F) {
  annot <- read.table(annotation_file, h=T,stringsAsFactors=F, sep = "\t")
  annot$ID <- paste(annot$Chr, annot$Start, annot$Ref, annot$Alt, sep="_")
  
  tmp <- annot[, c(1:56, ncol(annot))]
  annot <- unique(tmp)
  
  VAFs$ID <- rownames(VAFs)
  pass_annot <- merge(annot, VAFs[germline,], by="ID")
  write.csv(pass_annot, paste(output_dir, patient, "_passing_snps_annotated.csv", sep=""), quote = F, row.names = F)
  
}


## Next make a binary genotype matrix, both for passing and failing mutations
## Entries with 0.5 will not be used for treebuilding
genotype_bin <- as.matrix(passing_vafs)
genotype_bin[genotype_bin<gt_matrix_lower]=0
genotype_bin[genotype_bin>=gt_matrix_upper]=1
genotype_bin[genotype_bin>0&genotype_bin<1]=0.5

## Set mutations that shouldn't be counted to 'maybe'
genotype_bin[genotype_bin==1 & (!depth[germline,] | !mutant_reads[germline,])] <- 0

write.table(genotype_bin, paste(output_dir, patient, "_passing_snps_and_indels_genotype_binary_matrix.txt", sep=""), quote = F, row.names = T)

genotype_bin_snps <- subset(genotype_bin, rownames(genotype_bin) %in% rownames(NR_snp))
genotype_bin_indel <- subset(genotype_bin, rownames(genotype_bin) %in% rownames(NR_indels))
write.table(genotype_bin_snps, paste(output_dir, patient, "_passing_snps_genotype_binary_matrix.txt", sep=""), quote = F, row.names = T)
write.table(genotype_bin_indel, paste(output_dir, patient, "_passing_indels_genotype_binary_matrix.txt", sep=""), quote = F, row.names = T)

if(F) {
  genotype_bin_fail <- as.matrix(failing_vafs)
  genotype_bin_fail[genotype_bin_fail<gt_matrix_lower]=0
  genotype_bin_fail[genotype_bin_fail>=gt_matrix_upper]=1
  genotype_bin_fail[genotype_bin_fail>0&genotype_bin_fail<1]=0.5
  
  write.table(genotype_bin, paste(output_dir, patient, "_failed_snps_and_indels_genotype_binary_matrix.txt", sep=""), quote = F, row.names = T)
  
  
  png(paste(output_dir, patient, "_heatmap_passing_snps_and_indels.png", sep=""))
  heatmap(genotype_bin,scale='none',
          col=c("aliceblue","lightblue","steelblue"),mar=c(8,8))
  dev.off()
  
  png(paste(output_dir, patient, "_heatmap_passing_indels.png", sep=""))
  heatmap(genotype_bin_indel,scale='none',
          col=c("aliceblue","lightblue","steelblue"),mar=c(8,8))
  dev.off()
  
  #png(paste(output_dir, patient, "_heatmap_failed_snps_and_indels.png", sep=""))
  #heatmap(genotype_bin_fail,scale='none', col=c("aliceblue","lightblue","steelblue"),mar=c(8,8))
  #dev.off()
}


