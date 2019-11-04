#!/bin/bash

## Usage: This script is for making a consensus phylogenic tree for all the crypts of an individual. All the steps in
##        the script can be run individually, but this brings is all together for convenience.

## Before making a tree, here's what you need to do:
## 1. Run cgpVAF for all mutations passing Mathij's lmc-filters.
## 2. Run beta-binomial filters

## The steps in the tree making process are:
## 1. Read in the binary genotype matrix from beta-binomial filters script.
## 2. Remove crypts which are duplicates or should be excluded for any reason. Make the tree with and without.
## 3. Make a fasta file with all the sequences.
## 4. Run MPBoot. Do 3 and 4 separately for subs and indels.
## 5. Assign mutations to the tree branches. Write out the mutations for each branch for signature extraction.
## 6. Adjust the branch lengths

## Post Signature extraction - will be a separate script:
## 7. Plot trees with signatures
## 8. Make trees with branch lengths scaled to adjusted values of signature 1.

## Input: I expect the pileups to have been run and I expect the binomial_dir to contain a number of files...

## for i in 29 33 36 46 48 52; do bsub -o logs/patient${i}_indel.out -e logs/patient${i}_indel.err -M1000 -R'span[hosts=1] select[mem>1000] rusage[mem=1000]' "bash ~/phd/somatic_ibd_p1/phylogenics/make_MPBoot_tree_master.sh patient${i} indel"; done

pileup_dir=/lustre/scratch114/projects/crohns/somatic_ibd_p1/pileups/
binomial_dir=/lustre/scratch114/projects/crohns/somatic_ibd_p1/qc/binomial_filters/
script_dir=/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/phylogenics/

patient=$1
mutType=$2        ## snv or indel
MPBoot_dir=/lustre/scratch114/projects/crohns/somatic_ibd_p1/phylogenics/MPBoot/${mutType}/


mkdir -p ${MPBoot_dir}consensus_trees
mkdir -p ${MPBoot_dir}input_fasta_files
mkdir -p ${MPBoot_dir}binary_gt_matrices
mkdir -p ${MPBoot_dir}branch_mut_assignment
mkdir -p ${MPBoot_dir}tree_plots
mkdir -p ${MPBoot_dir}adj_mutCounts


if [ "$mutType" == "snv" ]; then
    echo "Making fasta file:"
    echo "/software/R-3.4.0/bin/Rscript ${script_dir}make_fasta_and_binary_gt_matrix.r ${patient} ${binomial_dir}${patient}/${patient}_passing_snps_genotype_binary_matrix.txt ${MPBoot_dir}input_fasta_files/"
    /software/R-3.4.0/bin/Rscript ${script_dir}make_fasta_and_binary_gt_matrix.r ${patient} ${mutType} \
    ${binomial_dir}${patient}/${patient}_passing_snps_genotype_binary_matrix.txt ${MPBoot_dir}input_fasta_files/

    echo "Running MPBoot"
    echo "${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}"
    ${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}consensus_trees/ ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}

    echo "Processing trees:"
    echo "/software/R-3.4.0/bin/Rscript ${script_dir}process_trees.R ${patient} ${binomial_dir}${patient}/${patient}_passing_snps_genotype_binary_matrix.txt ${MPBoot_dir}consensus_trees/${patient}/${patient}.treefile ${binomial_dir}${patient}/${patient}_NR_pass.txt ${binomial_dir}${patient}/${patient}_NV_pass.txt ${binomial_dir}${patient}/${patient}_medVAF_and_Unadj_MutCount.txt ${MPBoot_dir}"

    /software/R-3.4.0/bin/Rscript ${script_dir}process_trees.R ${patient} ${binomial_dir}${patient}/${patient}_passing_snps_genotype_binary_matrix.txt \
    ${MPBoot_dir}consensus_trees/${patient}/${patient}.treefile ${binomial_dir}${patient}/${patient}_NR_pass.txt \
    ${binomial_dir}${patient}/${patient}_NV_pass.txt ${binomial_dir}${patient}/${patient}_medVAF_and_Unadj_MutCount.txt \
    ${MPBoot_dir} ${mutType}
fi

if [ "$mutType" == "indel" ]; then
    echo "Processing indels for patient ${patient}"

    echo "Making fasta file:"
    /software/R-3.4.0/bin/Rscript ${script_dir}make_fasta_and_binary_gt_matrix.r ${patient} ${mutType} \
    ${binomial_dir}${patient}/${patient}_passing_indels_genotype_binary_matrix.txt ${MPBoot_dir}input_fasta_files/

    echo "Running MPBoot"
    echo "${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}"
    ${script_dir}run_MPBoot.sh ${MPBoot_dir} ${MPBoot_dir}consensus_trees/ ${MPBoot_dir}input_fasta_files/${patient}_trees.fa ${patient}

    echo "Processing trees:"
    /software/R-3.4.0/bin/Rscript ${script_dir}process_trees.R ${patient} ${binomial_dir}${patient}/${patient}_passing_indels_genotype_binary_matrix.txt \
    ${MPBoot_dir}consensus_trees/${patient}/${patient}.treefile ${binomial_dir}${patient}/${patient}_NR_pass.txt \
    ${binomial_dir}${patient}/${patient}_NV_pass.txt ${binomial_dir}${patient}/${patient}_medVAF_and_Unadj_MutCount.txt \
    ${MPBoot_dir} ${mutType}
fi

exit $?