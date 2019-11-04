
## The purpose of this script is to create a fasta file from all the mutations found for a patient that can then be
## used as input for MPBoot to generate a consensus phylogenic tree for the patient. 
## Based heavily on a script from Tim Coorens


options(stringsAsFactors = F)

if(file.exists("/nfs/users/nfs_s/so11/phd/doIexist")) {
  ## Script is running on the farm
  lib_location="/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm3/"
  require(ape, lib.loc=lib_location)
  require(seqinr, lib.loc=lib_location)
  .libPaths( c(lib_location, .libPaths()) )
} else {
  ## Running the script locally - Testing
  library(ape)
  library(seqinr)
}


args = commandArgs(TRUE)
patientID = args[1]
mutType=args[2]
binary_matrix = args[3]
fasta_output_dir= args[4]

#patientID <- "patient30"
#binary_matrix <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/qc/binomial_filters/patient30/patient30_passing_snps_and_indels_genotype_binary_matrix.txt"


genotype_bin <- read.table(binary_matrix, h=T)

# Create a fasta sequence to use as input for MPBoot
dna_strings = list()
Muts=rownames(genotype_bin)
samples <- gsub("_VAF", "", colnames(genotype_bin))

## Split the string like this to accommodate indels.
x <- unlist(strsplit(Muts, split ="_"))
Ref <- x[seq(3, length(x), 4)]
Alt <- x[seq(4, length(x), 4)]

## A hacky way to make the fasta string for indels but it works.
if(mutType=="indel") {
  Ref = rep("A",nrow(genotype_bin))
  Alt = rep("T",nrow(genotype_bin))
}

dna_strings <- list()
dna_strings[1]=paste(Ref,sep="",collapse="")
for (k in 1:length(samples)){
  Mutations = Ref
  Mutations[genotype_bin[,k]==0.5] = '?'
  Mutations[genotype_bin[,k]==1] = Alt[genotype_bin[,k]==1]
  dna_string = paste(Mutations,sep="",collapse="")
  dna_strings[k+1]=dna_string
}
names(dna_strings)=c("Ancestral",samples)


write.fasta(dna_strings, names=names(dna_strings), paste(fasta_output_dir, patientID, "_trees.fa", sep=""))




