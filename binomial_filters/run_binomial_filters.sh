#!/bin/bash

## The purpose of this script is to take in a crypt_metadata file and run the beta-binomial filters for all crypts
## in the file.
## The script will also plot heatmaps of the VAFs, print the global VAFs to a file

## 1. Merge the pileups for indels and snps and remove the vcf header.
## 2. Make beta-binomial script aware of patient sex. Make aware of the nr of crypts so it can skip if too few crypts.
## 3. Write out a table of global vafs for each mutation.
## 4. Write out A VAF heatmap for the patient and a VAF histogram
## 5. Write out a VAF heatmap of filtered mutations for the patient.
## 6. Make a file containing all mutations that pass for that patient.
## 7. Write a binary genotype matrix for the patient.

## Do I want to write a script to write on giant VAF histograms for all patients? Is this a way to visually compare
## clonality across IBD and control cohorts and possibly across disease subtypes and biopsy type.

## Make a script to filter mutations that occur in patients called using the same un-matched normal? Germline from P14
## called as somatic mutations?

## Note: This script only needs the massive amount of memory to plot the histograms in the end. If you don't care about that,
## the jobs can run with much less memory and be allowed to fail at the end.

## lcm_filtering_master.sh -m crypt_metaData -p patient_metaData [-s patient_subset=""] [-o output_dir=$PWD/] [-e clusterErrorOutputDir=output_dir/logs/] [-n nr_crypts_threshold=6] [-c cut_off]

## ~/phd/somatic_ibd_p1/qc/binomial_filters/run_binomial_filters.sh -m ~/phd/somatic_ibd_p1/sample_info/crypt_metadata_library_bank.txt -p ~/phd/somatic_ibd_p1/sample_info/patient_metaData.txt -s patient30


while [[ $# > 1 ]]
do
key="$1"

case $key in
    -m|--crypt_metaData)
    crypt_metaData="$2"
    shift
    ;;
    -p|--patient_metaData)
    patient_metaData="$2"
    shift
    ;;
    -o|--output|--output-dir)
    output_dir="$2"
    shift
    ;;
    -s|--patient_subset)
    patient_subset="$2"
    shift
    ;;
    -e|--clusterErrorOutput|--error_and_logfiles)
    clusterErrorOutput="$2"
    shift
    ;;
    -c|--cut_off)
    cut_off="$2"
    shift
    ;;
    -n|--nr_crypts_threshold)
    nr_crypts_threshold="$2"
    shift
    ;;
    --default)
    DEFAULT=YES
    shift
    ;;
    *)
    ;;
esac
shift
done

today=$( date | awk '{print $1"_"$2"_"$3}' )

echo "Binomial filtering of somatic mutations. log-file:" > logFile_qc_${today}.txt
date >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt
echo "Running ${worker_script}" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt

## Check input parameters and set defaults:
################################################
echo "Checking all input parameters and setting default values. Writing to logFile_qc_${today}.txt"
echo "Checking all input parameters and setting default values:" >> logFile_qc_${today}.txt
echo "---------------------------------------------------------" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt


if [ -z "$crypt_metaData" ]; then
    echo "Error: No crypt_metaData file specified. Please give the path to a file on the format described in the help message (Run this script with -h)" >> logFile_qc_${today}.txt
    echo "Error: No crypt_metaData file specified. Quitting" >> logFile_qc_${today}.txt
    exit 1
fi

if [ -z "$patient_metaData" ]; then
    echo "Error: No patient_metaData file specified. Please give the path to a file on the format described in the help message (Run this script with -h)" >> logFile_qc_${today}.txt
    echo "Error: No patient_metaData file specified. Quitting" >> logFile_qc_${today}.txt
    exit 1
fi




if [ ! -d "$output_dir" ]; then
    if [[ "$PWD" = */lustre/* ]]; then
        echo "No output directory specified or output directory does not exist. Writing output to the current working directory ${PWD}." >> logFile_qc_${today}.txt
        echo "" >> logFile_qc_${today}.txt
        echo "" >> logFile_qc_${today}.txt
        output_dir=${PWD}/
    else
        echo "Error: No output directory specified and the current working directory is not on lustre. Because this script
        interacts with the cluster, please use the -o option to specify an output directory on lustre." >> logFile_qc_${today}.txt
        echo "Quitting" >> logFile_qc_${today}.txt
        exit 1
    fi
else
    if [[ "$output_dir" = */lustre/* ]]; then
        echo "Writing output to ${output_dir}" >> logFile_qc_${today}.txt
        echo "" >> logFile_qc_${today}.txt
        echo "" >> logFile_qc_${today}.txt
    else
        echo "The specified output directory, ${output_dir}, is not on lustre. Because this script interacts with the
        cluster, please specify an output directory on lustre" >> logFile_qc_${today}.txt
        echo "Quitting" >> logFile_qc_${today}.txt
        exit 1
    fi
fi


if [ -z "$clusterErrorOutput" ]; then
    echo "No clusterErrorOutput directory specified. Writing log files under ${output_dir}logs/ by default." >> logFile_qc_${today}.txt
    echo "" >> logFile_qc_${today}.txt
    echo "" >> logFile_qc_${today}.txt
    mkdir -p ${output_dir}logs/
    clusterErrorOutput=${output_dir}logs/
fi
if [ ! -d "$clusterErrorOutput" ]; then
    echo "Error: The specified clusterErrorOutput directory, ${clusterErrorOutput} does not exist. "
    echo "" >> logFile_qc_${today}.txt
    echo "" >> logFile_qc_${today}.txt
    echo "Quitting" >> logFile_qc_${today}.txt
    exit 1
fi

if [ -z "$cut_off" ]; then
    cut_off=3
fi

if [[ ! $cut_off =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Invalid value for the cut_off parameter, ${cut_off}. Make sure this is a positive integer. The default value is 3." >> logFile_qc_${today}.txt
    echo "Quitting" >> logFile_qc_${today}.txt
    exit 1
fi

if [ -z "$nr_crypts_threshold" ]; then
    nr_crypts_threshold=5
fi

if [[ ! $nr_crypts_threshold =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Invalid value for the nr_crypts_threshold parameter, ${nr_crypts_threshold}. Make sure this is a positive integer. The default value is 6." >> logFile_qc_${today}.txt
    echo "Quitting" >> logFile_qc_${today}.txt
    exit 1
fi


## some hard-coded variables:
pileup_dir=/lustre/scratch114/projects/crohns/somatic_ibd_p1/pileups/
script_dir=/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/qc/binomial_filters/
lcm_filter_dir=/lustre/scratch114/projects/crohns/somatic_ibd_p1/qc/lcm_filter/


# Writing all input parameters to the log file
##############################################
echo "Parameters:" >> logFile_qc_${today}.txt
echo "-m crypt_meta data: $crypt_metaData" >> logFile_qc_${today}.txt
echo "-p patient_meta data: $crypt_metaData" >> logFile_qc_${today}.txt
echo "-o output directory: $output_dir" >> logFile_qc_${today}.txt
echo "-s patient_subset: ${patient_subset}" >> logFile_qc_${today}.txt
echo "-e clusterErrorOoutput: ${clusterErrorOutput}" >> logFile_qc_${today}.txt
echo "pileup_dir: ${pileup_dir}" >> logFile_qc_${today}.txt
echo "script_dir: ${script_dir}" >> logFile_qc_${today}.txt
echo "lcm_filter_dir: ${lcm_filter_dir}" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt
echo "" >> logFile_qc_${today}.txt

# Create a list of subjects
awk '$8=="WGS" && ($16=="NA" || $16=="Technical_Duplicate") {print $6}' < ${crypt_metaData} | sort -u > ${output_dir}patientList_binomialFiltering.tmp
if [ ! -z "$patient_subset" ]; then
    echo ${patient_subset} | awk -F, -v OFS="\n" '{$1=$1; print}' > ${output_dir}patient_subset.tmp
    grep -w -f ${output_dir}patient_subset.tmp ${output_dir}patientList_binomialFiltering.tmp > ${output_dir}tmp
    mv ${output_dir}tmp ${output_dir}patientList_binomialFiltering.tmp
fi


while read patient; do

    mkdir -p ${output_dir}${patient}
    # Is this IBD patient or control patient?
    cohort=$( awk -v PAT="$patient" '$6==PAT {print $18}' < ${crypt_metaData} | sort -u )

    ## What is the sex of the patient?
    sex=$( grep -w ${patient} ${patient_metaData} | awk '{print $3}' )
    if [ "$sex" == "M" ]; then
        sex=male
    else
        sex=female
    fi

    ## Remove the header lines and the comment char from both snp and indel input files
    cat < ${pileup_dir}${cohort}/${patient}/snp/*_snp_vaf.tsv | tail -n+2 | grep -v -e '^\##' | sed 's/\#//g' > ${pileup_dir}${cohort}/${patient}/snp/${patient}_snp_binomial_input.txt
    cat < ${pileup_dir}${cohort}/${patient}/indel/*_indel_vaf.tsv | grep -v -e '^\##' | sed 's/\#//g' > ${pileup_dir}${cohort}/${patient}/indel/${patient}_indel_binomial_input.txt
    ## Add this to deal with the homopolymers
    #cat ${pileup_dir}${cohort}/${patient}/indel_wo_homopolymers/*indel_vaf.tsv | grep -v -e '^\#' >> ${pileup_dir}${cohort}/${patient}/indel/${patient}_indel_binomial_input.txt

    ## Combine the annotation files from Mathij's filters. I want to merge with these to retain the annotation.
    cat ${lcm_filter_dir}${cohort}/${patient}/*_complete_final_allinfo_3.txt | awk 'NR==1 {print} $1!="Chr" {print}' > ${lcm_filter_dir}${cohort}/${patient}/${patient}_allinfo_combined.txt

    ## Submit R-script to the cluster
    echo "" >> logFile_qc_${today}.txt
    echo "" >> logFile_qc_${today}.txt
    echo "Submitted the following to the cluster:" >> logFile_qc_${today}.txt
    echo "    bsub -o ${clusterErrorOutput}/${patient}_binomial_log -e ${clusterErrorOutput}/${patient}_binomial_err -R'select[mem>20000] rusage[mem=20000]' -M20000 \
    '/software/R-3.4.2/bin/Rscript ${script_dir}binomial_filters_worker.R ${pileup_dir}${cohort}/${patient}/snp/${patient}_snp_binomial_input.txt ${pileup_dir}${cohort}/${patient}/indel/${patient}_indel_binomial_input.txt \
    ${lcm_filter_dir}${cohort}/${patient}/${patient}_allinfo_combined.txt ${sex} ${cut_off} ${nr_crypts_threshold} ${script_dir} ${output_dir}${patient}/ ${patient}'"  >> logFile_qc_${today}.txt

    bsub -q normal -o ${clusterErrorOutput}/${patient}_binomial_log -e ${clusterErrorOutput}/${patient}_binomial_err -R"select[mem>20000] rusage[mem=20000]" -M20000 \
    "/software/R-3.4.2/bin/Rscript ${script_dir}binomial_filters_worker.R ${pileup_dir}${cohort}/${patient}/snp/${patient}_snp_binomial_input.txt ${pileup_dir}${cohort}/${patient}/indel/${patient}_indel_binomial_input.txt \
    ${lcm_filter_dir}${cohort}/${patient}/${patient}_allinfo_combined.txt ${sex} ${cut_off} ${nr_crypts_threshold} ${script_dir} ${output_dir}${patient}/ ${patient}"



done < ${output_dir}patientList_binomialFiltering.tmp



exit $?