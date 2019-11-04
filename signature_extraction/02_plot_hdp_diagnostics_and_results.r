

hdp_posterior_dir="/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/"
output_dir <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/"
# read in PCAWG sigs
all_sigs <- read.csv("/Users/so11/phd/so11_nfs/somatic_ibd_p1/resources/pcawk_signatures/sigProfiler_SBS_signatures_summing_to_one.csv", header = T, row.names = 1, stringsAsFactors = F)

library(hdp)

load(paste(output_dir, "branch_mutational_matrix.txt", sep=""))
tnt <- mut_mat_t
newColNames <- paste(substr(colnames(mut_mat_t), 3,3), ".", substr(colnames(mut_mat_t), 5,5), ".in.", substr(colnames(mut_mat_t), 1,1), substr(colnames(mut_mat_t), 3,3), substr(colnames(mut_mat_t), 7,7), sep="")
colnames(tnt) <- newColNames
patients <- sapply(strsplit(rownames(tnt), "_"), "[[", 1)



worked <- c(1:20)

for(i in worked) {
  load(paste(hdp_posterior_dir, "hdpout_", i, ".RData", sep=""))
}

chlist <- vector("list", length(worked))
for(i in worked) {
  chlist[[i]] <- get(paste("hdp_", i, sep=""))
}


chlist <- vector("list", length(worked)-1)
chlist[[1]] <- hdp_1
chlist[[2]] <- hdp_2
chlist[[3]] <- hdp_3
chlist[[4]] <- hdp_4
chlist[[5]] <- hdp_5
chlist[[6]] <- hdp_6
chlist[[7]] <- hdp_7
chlist[[8]] <- hdp_8
chlist[[9]] <- hdp_9
chlist[[10]] <- hdp_10
chlist[[11]] <- hdp_11
chlist[[12]] <- hdp_12
chlist[[13]] <- hdp_13
chlist[[14]] <- hdp_14
chlist[[15]] <- hdp_15
chlist[[16]] <- hdp_16
chlist[[17]] <- hdp_18
chlist[[18]] <- hdp_19
chlist[[19]] <- hdp_20


luad_multi <- hdp_multi_chain(chlist)

# plot diagnostics
par(mfrow=c(4,5), mar=c(5,4,1,1))
lapply(chains(luad_multi), plot_lik, bty='L', start=1000)
lapply(chains(luad_multi), plot_lik, bty='L')
lapply(chains(luad_multi), plot_numcluster, bty='L')
lapply(chains(luad_multi), plot_data_assigned, bty='L')



## extract consensus components / signatures
luad_multi <- hdp_extract_components(luad_multi, cos.merge = 0.85)
luad_multi
numcomp(luad_multi)
prop.ex(luad_multi)

par(mfrow=c(1,1), mar=c(4, 4, 2, 2))

# plot number of mutations per signature (one dot per posterior sample)
plot_comp_size(luad_multi, bty="L")

# plot components / signatures
# pick your colours
#mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')
## Use same colours as COSMIC
mut_colours <- c("#16bdebff", "#000000ff", "#e22926ff", "#999999ff", "#9fce67ff", "#ecc6c5ff")
# labels along bototm x-axis
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
# group labels along the top (and controls colour grouping)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

par(mfrow=c(4,4), mar=c(4, 4, 2, 2))
## 
plot_comp_distn(luad_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)

## Once you've figured out which signature is which then you can
## make the plot with signature names, but this will have to be done
## manually as the number and order of signatures will change.
if(FALSE) {
  plot_comp_distn(luad_multi, cat_names=trinuc_context,
                  grouping=group_factor, col=mut_colours,
                  col_nonsig="grey80", show_group_labels=TRUE,
                  plot_title = c("Uncertainty","Signature1","Signature5", "Chemotherapy signature", "SBSA - Discovered by Henry",
                                 "Signature18","SBSB (Henry)", "Purine signature", "Noise?", "SBSC (Henry)", 
                                 "Signature13", "New signature?", "Signature2",
                                 "Signature17b", "Signature17a", "Signature45"))
  
}

par(mfrow=c(1,1), mar=c(4, 4, 2, 2))
## The object contains 1 root node, one node that is the parent of all data, 1 pseuodo node for each prior
## and one node for each patient. 
## 1+1+length(gdsigs) + length(unique(patients))
## The first node containing sample information then is the next node after that. 
## This number will change as more samples and a different number of priors are used so beware. 
first_sample_node <- 1+1+length(gdsigs) + length(unique(patients)) + 1
last_sample_node <- length(comp_dp_distn(luad_multi)$cred.int)

## Need to come up with a better solution for the color palette. Need a palette that can handle more colours. 
plot_dp_comp_exposure(luad_multi, main_text="Signatures - All patients",
                      dpindices=first_sample_node:last_sample_node,
                      col=c(RColorBrewer::brewer.pal(12, "Set3"), "navy", "forestgreen", "red", "yellow"),
                      incl_nonsig=FALSE,
                      ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                      leg.title = 'Signature',incl_numdata_plot = T, dpnames = rownames(tnt),
                      las=2,mar=c(3,4,2,0.5))


## Plot the control and IBD cohorts separately
## Assuming here that there's a continous string of IBD patient samples - this may not be the case. 
first_ibd_node <- min(grep("patient", rownames(tnt)))

plot_dp_comp_exposure(luad_multi, main_text="Signatures - IBD patients",
                      dpindices=first_ibd_node:last_sample_node,
                      col=c(RColorBrewer::brewer.pal(12, "Set3"), "navy", "forestgreen", "red", "yellow"),
                      incl_nonsig=FALSE,
                      ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                      leg.title = 'Signature',incl_numdata_plot = T, dpnames = rownames(tnt)[(first_ibd_node:last_sample_node)-first_sample_node],
                      las=2,mar=c(3,4,2,0.5))

plot_dp_comp_exposure(luad_multi, main_text="Signatures - Controls",
                      dpindices=first_sample_node:(first_ibd_node-1),
                      col=c(RColorBrewer::brewer.pal(12, "Set3"), "navy", "forestgreen", "red", "yellow"),
                      incl_nonsig=FALSE,
                      ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                      leg.title = 'Signature',incl_numdata_plot = T, dpnames = rownames(tnt)[(first_sample_node:(first_ibd_node-1))-(first_sample_node-1)],
                      las=2,mar=c(3,4,2,0.5))

plot_patient_dp_comp <- function(hdp_multi= luad_multi, patient, first_sample_node, trinuc_matrix=tnt, 
                                 incl_nonSig=F, incl_numData=T) {
  firstNode=min(grep(paste(patient, "_", sep=""), rownames(trinuc_matrix))) + first_sample_node -1
  lastNode=max(grep(paste(patient, "_", sep=""), rownames(trinuc_matrix))) + first_sample_node -1
  
  plot <- plot_dp_comp_exposure(luad_multi, main_text=patient,
                                dpindices=firstNode:lastNode,
                                col=c(RColorBrewer::brewer.pal(12, "Set3"), "navy", "forestgreen", "red", "yellow"),
                                incl_nonsig=incl_nonSig,
                                ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                                leg.title = 'Signature',incl_numdata_plot = incl_numData, dpnames = rownames(tnt)[(firstNode:(lastNode))-(first_sample_node -1)],
                                las=2,mar=c(3,4,2,0.5))
  print(plot)
}



# to extract relative contributions
toutmeans <- (comp_dp_distn(luad_multi)$mean)
toutmeans <- toutmeans[first_sample_node:(nrow(toutmeans)),]
rownames(toutmeans) <- rownames(tnt)

write.table(data.frame(toutmeans),"/Users/so11/phd/so11_lustre/somatic_ibd_p1/signature_extraction/hdp/snvs/branch_exposures_w_prior.txt", row.names = T, quote=F)
