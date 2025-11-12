#!/bin/bash

# TO NOTE: had initially written this to make two filesets - one with
# genotypes and one that saved dosages. Dosages version was too large
# to be feasible. If return with time to troubelshoot, should revisit.

# TO NOTE: this script has one option that uses PLINK and another that
# does  not, so preserves all dosage information in imputed VCF (helpful
# for keeping HDS after Michigan HLA imputation). This also means that
# the non-PLINK option is much slower and should only be run if more
# than GT is needed.

set -e
set -u

while getopts n:r:m:p:d:o: opt; do
   case "${opt}" in
      n) chr=${OPTARG};;
      r) rsq=${OPTARG};;
      m) maf=${OPTARG};;
      p) plink_prefix=${OPTARG};;
      d) dir=${OPTARG};;
      o) option=${OPTARG};;
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done


# Make out dir
in_dir="${dir}/imputed"
out_dir="${dir}/imputed_clean_maf${maf}_rsq${rsq}"
mkdir -p "${dir}/imputed_clean_maf${maf}_rsq${rsq}"

# Set filter
to_filt="((INFO/TYPED = 1 | (INFO/IMPUTED = 1 & INFO/R2 > ${rsq})) & INFO/MAF > ${maf})"

# Filter INFO file, for smaller output with kept variables
    # ((typed OR (imputed & R2)) & MAF)
bcftools filter -i \
    "$to_filt" \
    "${in_dir}/chr${chr}.info.gz" -o "${out_dir}/chr${chr}_clean.info"

# Write out list of RSIDs want to keep
bcftools query -f '%ID\n' \
    "${out_dir}/chr${chr}_clean.info" > "${out_dir}/chr${chr}_maf${maf}_rsq${rsq}_snps.txt"

# If option set to use PLINK for only GT
if [ "$option" = "gt" ]; then
    # Filter VCF to these IDs using PLINK, keeping GT
    if [ "$chr" = "X" ]; then
        # If chrX, PLINK throws error without sex information
        plink2 --vcf "${in_dir}/chr${chr}.dose.vcf.gz" \
            --export vcf 'bgz' \
            --extract "${out_dir}/chr${chr}_maf${maf}_rsq${rsq}_snps.txt" \
            --id-delim '_' \
            --fam "${plink_prefix}.fam" \
            --out "${out_dir}/tmp_chr${chr}_clean"
    else
        plink2 --vcf "${in_dir}/chr${chr}.dose.vcf.gz" \
            --export vcf 'bgz' \
            --extract "${out_dir}/chr${chr}_maf${maf}_rsq${rsq}_snps.txt" \
            --out "${out_dir}/tmp_chr${chr}_clean"
    fi

# if want to keep all dosage information, use bcftools (much slower)
elif [ "$option" = "all" ]; then
    # Filter VCF to these IDs using bcftools
    # NOTE: will add sex information at conversion to PLINK in later script
    bcftools view --include ID==@"$snp_list" "${in_dir}/chr${chr}.dose.vcf.gz" \
        -Oz -o "${out_dir}/tmp_chr${chr}_clean.vcf.gz"
fi

# Finally, clean up RSIDs that may have appeared more than once, mainly '.' IDs.
bcftools filter -i \
    "$to_filt" \
    "${out_dir}/tmp_chr${chr}_clean.vcf.gz" -o "${out_dir}/chr${chr}_clean.vcf.gz"
tabix "${out_dir}/chr${chr}_clean.vcf.gz"

# Clean up
rm -f ${out_dir}/tmp_chr${chr}_*
