# Somatic evolution in non-neoplastic IBD-affected colon. Olafsson et al.

## A guide to reproducing the results of the manuscript and using the data

### Overview of resources:

Scripts required to reproduce the analyses of the manuscript are available at:
https://github.com/Solafsson/somaticIBD
See the Key Resources table in the STAR methods for a complete list of all software and from where it can be downloaded.

Raw sequencing data can be accessed from the European Genome-phenome Archive (EGA) with the accession code EGAD00001006061

The raw data for the control cohort can be accessed from the European Genome-phenome Archive (EGA) with accession codes EGAD00001004192 and EGAD00001004193.

Mendeley Data repository: http://dx.doi.org/10.17632/x3vsxpspn4.2

I have deposited several intermediate files in Mendeley Data. These include:
* All_IBD_cohort_mutations_mapped_to_branches.txt and the corresponding file for controls contain the full lists of mutations mapped to branches of the phylogenetic trees. The sampleID columns refer to node-numbers of the tree for each particular patient.
* Phylogenetic_trees_SBS/ contains the tree files generated by MPBoot. The length of the branches represents the number of substitutions on the branch. I only used the substitutions to create the trees and then mapped the indels to the SBS-trees. The node-numbers in these trees are the same as the SampleIDs in All_IBD_cohort_mutations_mapped_to_branches.txt. I have provided both the trees with the unadjusted branch lengths and the trees with branch lengths adjusted for coverage and VAF as described in the methods section of the paper.
* Pile_ups_read_counts/ contains for each crypt from every patient the total number of reads covering any mutation detected in crypts from that individual, and a second matrix containing the number of reads supporting the ALT allele.
* Signature_extraction/ contains the data needed to reproduce the signature extraction. For snvs and indels, it contains the PCAWG signatures used as priors and the 96 or 83 class mutation matrices that serve as input. Signature extraction was performed on a per-branch basis and so each line of the input refers to a branch of the phylogenetic trees in Phylogenetic_trees_SBS/. I also give the output of the signature extraction on a per-branch basis.
* Spatial_relationship_matrices/ contains the physical distances between crypts in the unit of crypt-widths. Refer to the Histology_images/


### Mutation calling
If you don’t want to use the filtered mutation calls provided in the Mendeley directory and the supplementary tables you can follow the links given in the Key Resources table in the STAR methods to set up CaveMan, Pindel, ASCAT, BRASS and TraFiC-mem.
These are all part of the Sanger pipeline and so I can’t provide any advice on set-up or installation.

Some patients have a matched normal sample and some don’t. When a matched normal sample is available, run all algorithms using that sample as matched normal. Otherwise, run unmatched with P14B2_No2b as normal.

For ASCAT, BRASS and TraFiC-mem, you want to use a sample from the same individual as normal to remove germline variants. I generally picked a clonally unrelated sample with the highest coverage.


### Filtering
Mathijs wrote some custom filters both for CaveMan/Pindel and BRASS calls that are specially tailored to capture artefacts associated with the small input pipeline of the laser capture microscopy pipeline.
The BRASS filters are available here: https://github.com/MathijsSanders/AnnotateBRASS and the filters applied to SNVs and Indels are here: https://github.com/MathijsSanders/SangerLCMFiltering

After applying these basic filters, a large number of germline variants will remain. Use the scripts in the binomial_filters/ directory on GitHub to apply the binomial and beta-binomial filters described in the methods section. This will only remove a small number of variants from samples that were called using a matched normal but will be vital to remove germline variants from samples that underwent unmatched calling.

### Tree building
The next step is the tree building. Download MPBoot using the link in the Key Resources table and use the scripts in phylogenics/ in the github directory to build phylogenetic trees for all samples, to adjust the mutation burden of each crypt for coverage and VAF and to assign mutations to branches.

### Signature extraction
Signatures are extracted treating each branch with >50 mutations as a sample. Use the scripts provided in https://github.com/Solafsson/somaticIBD

### Selection analyses

I used the dndscv software (see https://github.com/im3sanger/dndscv). Please see the relevant script in https://github.com/Solafsson/somaticIBD


### Burden analyses

This includes regressing substitution and indel burdens against age and disease duration. Also includes regressions of signature burden with the same. Also includes association testing of SV status and driver carrier status with mutation burden. Please see the relevant script in https://github.com/Solafsson/somaticIBD
