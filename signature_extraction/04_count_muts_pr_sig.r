

library(ape)
library(ggtree)

patient_meta <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/sample_info/patient_metaData.txt", h=T)
crypt_meta <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt",h=T)
tree_dir="/Users/so11/phd/so11_lustre/somatic_ibd_p1/phylogenics/MPBoot/snv/consensus_trees/"
toutmeans <- read.table("/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/branch_exposures_w_prior.txt", h=T)

crypt_meta$patient_ID <- as.character(crypt_meta$patient_ID)
patients <- unique(crypt_meta$patient_ID[is.na(crypt_meta$excl_criteria) & crypt_meta$WGS_status=="WGS"])
## Some patients to exclude. O337 is cancer patient, patient7 is the same as patient4.
patients <- patients[!(patients %in% c("O337", "patient7", "patient30"))]

## Hard-coded. May need to change if you re-run the HDP
sigNames <- c("Unassigned", "SBS1","SBS5", "SBSA", "SBS18","SBS32", "SBSB","SBS1_noise","SBS17b", "SBS35", "SBS-N1",
              "SBS-N2", "SBS-N3", "SBS13", "SBS2", "SBS17a")

all_sig_muts <- data.frame()
for(patient in patients) {

  tree <- read.tree(paste(tree_dir, patient, "/", patient, "_MutAssign_truncBinom_adj.tree", sep=""))
  tree_original <- tree
  tree_df <- fortify(tree_original)

  ## Extract the branches for this patient. 
  b <- toutmeans[grep(paste(patient, "_", sep=""), rownames(toutmeans)),]
  b$node <- unlist(strsplit(rownames(b), split="_"))[c(FALSE, TRUE)]

  test <- merge(tree_df, b, by="node", all.x=T)
  
  ## For each signature, multiply the branch lengths with the fraction of mutations attributed to that signature
  ## in that branch. Then sum over all 'anchestor branches' to obtain total counts for each crypt. 
  for(i in 1:length(sigNames)) {
    tree_df <- fortify(tree_original)
    adj.branch <- test$branch.length

    adj.branch[!is.na(test[,(9+i)])] <- test$branch.length[!is.na(test[,(9+i)])]*test[!is.na(test[,(9+i)]), (9+i)]
    tree$edge.length=adj.branch[tree$edge[,2]]
    tree_df <- fortify(tree)
    if(i==1) {
      test_sigs <- tree_df[tree_df$isTip, c("label", "x")]
      
      ## Some branches have fewer than 50 mutations and so were not included in the signature extraction. 
      ## I want to find the total number of mutations belonging to such branches that precede each crypt
      ## and subtract that number from the number of mutations assigned to each signature. 
      no_extraction_branches <- adj.branch
      no_extraction_branches[!is.na(test[,(9+i)])] <- 0
      tree_tmp <- tree_original
      tree_tmp$edge.length=no_extraction_branches[tree_tmp$edge[,2]]
      tree_tmp_df <- fortify(tree_tmp)
      no_ext <- tree_tmp_df[tree_tmp_df$isTip, "x"]
    } else {
      test_sigs <- cbind(test_sigs, (tree_df[tree_df$isTip, c("x")]-no_ext))
    }
    #write.tree(tree, paste(tree_dir, patient, "/", patient, "_", sigNames[i], "_truncBinom_adj.tree", sep=""))
  }
  colnames(test_sigs) <- c("crypt_ID", sigNames)
  all_sig_muts <- rbind(all_sig_muts, test_sigs)

}


crypt_meta <- merge(crypt_meta, patient_meta, by="patient_ID")
crypt_meta <- merge(crypt_meta, all_sig_muts, by="crypt_ID", all.x=T)
crypt_meta <- subset(crypt_meta, WGS_status=="WGS" & type!="Matched_Normal" & !is.na("excl_criteria"))


## Plot a stacke barplot with signatures for UC and CD separately. 
####################################################################
library(reshape2)
library(ggplot2)
library(cowplot)

## Clean up the signatures
## There are some components from the HDP which we don't believe are real signatures. Lump these together with the
## Unassigned fraction. 
crypt_meta$Unassigned <- crypt_meta$Unassigned + crypt_meta$`SBS-N1` + crypt_meta$`SBS-N2` + crypt_meta$`SBS-N3`
crypt_meta$`SBS-N1`<- NULL
crypt_meta$`SBS-N2`<- NULL
crypt_meta$`SBS-N3`<- NULL

## HDP had problems deconvoluting some highly correlated signatures. Estimated SBS5 and SBS18 within SBS1
## in EM_to_share_out_sigs.r and now want to account for this. 
sbs1_decon <- read.table("/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/EM_reconstruction/hdp_SBS_broken_down_into_pcawg_1-5-18.txt", h=T)
crypt_meta$SBS5 <- crypt_meta$SBS5 + sbs1_decon[2,2]*crypt_meta$SBS1
crypt_meta$SBS18 <- crypt_meta$SBS18 + sbs1_decon[3,2]*crypt_meta$SBS1
crypt_meta$SBS1 <- sbs1_decon[1,2]*crypt_meta$SBS1

## Add N3 to the SBS1 component. We attribute this to inter-individual variance in SBS1. The second CpG column is 
## sometimes as high as the first one... 
crypt_meta$SBS1 <- crypt_meta$SBS1 + crypt_meta$SBS1_noise
crypt_meta$SBS1_noise <- NULL

sig_colours <- c("#F5BB00","#53ADFC","#FCECC9","#3BB273", "#DD3B46","#C46BAE", "#0566C6","#B5ABAB", 
                 "#F17300", "#A7CECB", "#2D2A32",
                 "#247BA0")

CD_crypts <- subset(crypt_meta, Disease=="CD")
CD_sigs <- CD_crypts[, c(1,2, (ncol(crypt_meta)-11):ncol(crypt_meta))]
CD_sigs <- CD_sigs[complete.cases(CD_sigs),]
CD_sigs$totalMuts <- apply(CD_sigs[,c(3:14)], 1, sum)
CD_sigs <- CD_sigs[order(CD_sigs$totalMuts), ]
CD_sigs$nr <- c(1:nrow(CD_sigs))

UC_crypts <- subset(crypt_meta, Disease=="UC")
UC_sigs <- UC_crypts[, c(1,2, (ncol(crypt_meta)-11):ncol(crypt_meta))]
UC_sigs <- UC_sigs[complete.cases(UC_sigs),]
UC_sigs$totalMuts <- apply(UC_sigs[,c(3:14)], 1, sum)
UC_sigs <- UC_sigs[order(UC_sigs$totalMuts), ]
UC_sigs$nr <- c(1:nrow(UC_sigs))

sig_CD_m <- melt(CD_sigs, id.vars=c("crypt_ID", "nr", "totalMuts", "patient_ID"))
sig_UC_m <- melt(UC_sigs, id.vars=c("crypt_ID", "nr", "totalMuts", "patient_ID"))
sig_CD_m$Disease <- "Crohn's disease"
sig_UC_m$Disease <- "Ulcerative colitis"
sigs_m <- rbind(sig_CD_m, sig_UC_m)

## Order patients by highest mutation burden:
sig_CD_m$patient_ID <- factor(sig_CD_m$patient_ID, levels=names(tapply(sig_CD_m$totalMuts, sig_CD_m$patient_ID, max))
                              [order(tapply(sig_CD_m$totalMuts, sig_CD_m$patient_ID, max), decreasing = T)])

sig_UC_m$patient_ID <- factor(sig_UC_m$patient_ID, levels=names(tapply(sig_UC_m$totalMuts, sig_UC_m$patient_ID, max))
                              [order(tapply(sig_UC_m$totalMuts, sig_UC_m$patient_ID, max), decreasing = T)])

## Add a little explanation to the signature names 
sig_CD_m$variable <- factor(sig_CD_m$variable, labels=c("Unassigned", "SBS1 (CpG deamin)","SBS5 (Clock-like)", "SBSA (Unknown)", "SBS18 (ROS)","SBS32 (Purine)", "SBSB (Unknown)", "SBS17b (Unknown)", "SBS35 (Platinum)", 
                                                                     "SBS13 (APOBEC)", "SBS2 (APOBEC)", "SBS17a (Unknown)"))

p1 <- ggplot(sig_CD_m, aes(fill=variable, y=value, x=reorder(crypt_ID, -nr))) + theme_bw() + 
  geom_bar(position="fill", stat="identity", width = 1) + theme(axis.text.x = element_blank()) + 
  facet_grid(~patient_ID, scale="free",space="free_x") + labs(x="Individual crypts from Crohn's disease", y="Proportion", fill="Mutation Signature") +
  theme(plot.margin=unit(c(-0.35,1,1,1), "cm")) + scale_fill_manual(values=sig_colours) + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())  +
  theme(axis.text = element_text(size=30), axis.title.y=element_text(size=30), axis.title.x=element_text(size=45)) + 
  theme(legend.text = element_text(size=30), legend.key.size = unit(1.5, "cm"), legend.title = element_text(size=36))



p2 <- ggplot(sig_CD_m, aes(x=reorder(crypt_ID, -nr), y =totalMuts/(length(sigNames)))) + theme_bw() + 
  geom_bar(stat="identity", width=1) + theme(axis.text.x = element_blank()) + 
  facet_grid(~patient_ID, scale="free",space="free_x") + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank()) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),plot.margin=unit(c(1,1,0.5,1), "cm")) + 
  labs(x="", y="Total subs") + theme(axis.text = element_text(size=30), axis.title.y=element_text(size=30)) 


p3 <- ggplot(sig_UC_m, aes(fill=variable, y=value, x=reorder(crypt_ID, -nr))) + theme_bw() + 
  geom_bar(position="fill", stat="identity", width = 1) + theme(axis.text.x = element_blank()) + 
  facet_grid(~patient_ID, scale="free",space="free_x") + labs(x="Individual crypts from Ulcerative colitis", y="Proportion", fill="Mutation Signature") +
  theme(plot.margin=unit(c(-0.35,1,1,1), "cm")) + scale_fill_manual(values=sig_colours) + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank())  +
  theme(axis.text = element_text(size=30), axis.title.y=element_text(size=30), axis.title.x=element_text(size=45)) + 
  theme(legend.position = "none")



p4 <- ggplot(sig_UC_m, aes(x=reorder(crypt_ID, -nr), y =totalMuts/(length(sigNames)))) + theme_bw() + 
  geom_bar(stat="identity", width=1) + theme(axis.text.x = element_blank()) + 
  facet_grid(~patient_ID, scale="free",space="free_x") + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank()) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),plot.margin=unit(c(1,1,0.5,1), "cm")) + 
  labs(x="", y="Total subs") + theme(axis.text = element_text(size=30), axis.title.y=element_text(size=30)) 



png("/Users/so11/phd/somatic_ibd_p1/manuscript/Figures/Figure3A_byPat.png", width = 2500, height = 1400 )
p <- plot_grid(p2, p1,p4,p3,  rel_heights = c(1,2), nrow=4, align = "v", axis="r")
print(p)
dev.off()


## Want to write out the number of mutations attributed to each signature.
## I will test these for association with disease duration in the lmm step.
SignatureExposure <- crypt_meta[is.na(crypt_meta$excl_criteria) & crypt_meta$WGS_status=="WGS" & crypt_meta$type!="Matched_Normal",
                                c(1, (ncol(crypt_meta)-11):ncol(crypt_meta))]

write.table(SignatureExposure, "/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/crypt_exposure_mutCounts_post_EM.txt", row.names = F, quote = F)





