#!/bin/bash

set -e
set -u

while getopts p:d: opt; do
   case "${opt}" in
      p) plink_prefix=${OPTARG};;
      d) dir=${OPTARG};;
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done


# Make file_list for input into bcftools concat
ls ${dir}/*_clean.vcf.gz | sort -V > "${dir}/imputed_cleaned_file_list.txt"
 
# Concatenate VCFs
bcftools concat \
   --file-list "${dir}/imputed_cleaned_file_list.txt" \
   --output-type z \
   --output "${dir}/chr_all_concat.vcf.gz"

# Convert to PLINK, adding sex info from pre_qc files
plink2 --vcf "${dir}/chr_all_concat.vcf.gz" \
   --make-pgen --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 10000 \
   --id-delim '_' \
   --update-sex "${plink_prefix}.fam" col-num=5 \
   --out "${dir}/chr_all_concat"
