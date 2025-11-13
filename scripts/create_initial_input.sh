#!/bin/bash

set -e
set -u

# Set defaults for optional args
id_keep_list=""
hwe_filt_list=""

while getopts p:o:c:n:b:t:s:k:h:l: opt; do
   case "${opt}" in
      p) plink_prefix=${OPTARG};;
      o) out_dir=${OPTARG};;
      c) code_dir=${OPTARG};;
      n) chr=${OPTARG};;
      b) orig_build=${OPTARG};;
      t) to_build=${OPTARG};;
      s) sexcheck=${OPTARG};;
      k) id_keep_list=${OPTARG};;
      h) hwe_filt_list=${OPTARG} ;;
      l) log_file=${OPTARG};;
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done

# Current code ALWAYS applies less strict HWE for chr6.
# Note that PLINK1.9 input is required.
# Note that id_keep_list and hwe_filt_list are optional.

#TODO: should this be split into different bash scripts? Different
# snakefile rules?
#TODO: better way to handle all the file copying here for different
# branching logic?

#### Pre-filtering ------------------------------------------------------------

# Function to check file exists, use for awk and other commands that could
# fail silently
check_file_exists() {
   if [ ! -f "$1" ]; then
      echo "Error: Required file '$1' not found." >&2
      exit 1
   fi
}

# Convert to array
chr=($chr)
sexcheck=($sexcheck)

# Set up log file
echo -e "Category\tDescription\tSamples\tSNPs" > "$log_file"
check_file_exists "${plink_prefix}.bim"
ct_var_plink_prefix=$(awk 'END{print NR}' "${plink_prefix}.bim")
ct_samp_plink_prefix=$(awk 'END{print NR}' "${plink_prefix}.fam")
echo -e "Summary\tRaw data counts\t${ct_samp_plink_prefix}\t${ct_var_plink_prefix}" >> "$log_file"

# Subset to list of IDs, if provided
plink_prefix_name=$(basename "${plink_prefix}")
if [ -n "$id_keep_list" ]; then
   id_keep_list_name=$(basename "$id_keep_list")
   echo "Optional: Filtering to IDs in given list: ${id_keep_list_name}."
   plink2 --bfile "$plink_prefix" \
      --keep "$id_keep_list" \
      --make-bed \
      --out "${out_dir}/${plink_prefix_name}_subset"

   check_file_exists "${out_dir}/${plink_prefix_name}_subset.fam"
   ct_subset=$(awk 'END{print NR}' "${out_dir}/${plink_prefix_name}_subset.fam")
   subset_rm=$((ct_samp_plink_prefix - ct_subset))

   echo -e "Pre-Filtering\tSubset based on given ID list\t(${subset_rm})\t()" >> "$log_file"

   # Copy to tmp file so can use same name if subset not done
   cp "${out_dir}/${plink_prefix_name}_subset.bim" "${out_dir}/tmp_subset.bim"
   cp "${out_dir}/${plink_prefix_name}_subset.bed" "${out_dir}/tmp_subset.bed"
   cp "${out_dir}/${plink_prefix_name}_subset.fam" "${out_dir}/tmp_subset.fam"
else
   plink2 --bfile $plink_prefix \
      --make-bed --out ${out_dir}/tmp_subset
   cp "${plink_prefix}.bim" "${out_dir}/tmp_subset.bim"
   cp "${plink_prefix}.bed" "${out_dir}/tmp_subset.bed"
   cp "${plink_prefix}.fam" "${out_dir}/tmp_subset.fam"
fi

# Remove SNPs missingness > 20%
plink2 --bfile "${out_dir}/tmp_subset" \
   --geno 0.2 \
   --make-bed \
   --out "${out_dir}/tmp_miss_rm"
check_file_exists "${out_dir}/tmp_subset.bim"
check_file_exists "${out_dir}/tmp_miss_rm.bim"
ct_orig=$(awk 'END{print NR}' "${out_dir}/tmp_subset.bim")
ct_miss=$(awk 'END{print NR}' "${out_dir}/tmp_miss_rm.bim")
miss_rm=$((ct_orig - ct_miss))
echo -e "Pre-Filtering\tRemove SNPs with missingness > 20%\t()\t($miss_rm)" >> "$log_file"

# Remove monomorphic SNPs
plink2 --bfile "${out_dir}/tmp_miss_rm" \
   --freq \
   --nonfounders \
   --out "${out_dir}/tmp_mono"

Rscript ${code_dir}/scripts/filt_mono.R \
   -a "${out_dir}/tmp_mono.afreq" \
   -l "$log_file"

plink2 \
   --bfile "${out_dir}/tmp_miss_rm" \
   --exclude "${out_dir}/tmp_mono_rm.txt" \
   --make-bed \
   --out "${out_dir}/tmp_mono_rm"

# Remove samples with call rate < 95%
plink2 --bfile "${out_dir}/tmp_mono_rm" \
   --mind 0.05 \
   --nonfounders \
   --make-bed \
   --out "${out_dir}/tmp_samp_rm"
check_file_exists "${out_dir}/tmp_mono_rm.fam"
check_file_exists "${out_dir}/tmp_samp_rm.fam"
ct_orig=$(awk 'END{print NR}' "${out_dir}/tmp_mono_rm.fam")
ct_samp=$(awk 'END{print NR}' "${out_dir}/tmp_samp_rm.fam")
samp_rm=$((ct_orig - ct_samp))
echo -e "Pre-Filtering\tRemove samples with missingness > 5%\t($samp_rm)\t()" >> "$log_file"

# Remove SNPs missingness > 5%
plink2 --bfile "${out_dir}/tmp_samp_rm" \
   --geno 0.05 \
   --nonfounders \
   --make-bed \
   --out "${out_dir}/tmp_miss_rm_2"
check_file_exists "${out_dir}/tmp_samp_rm.bim"
check_file_exists "${out_dir}/tmp_miss_rm_2.bim"
ct_orig=$(awk 'END{print NR}' "${out_dir}/tmp_samp_rm.bim")
ct_miss=$(awk 'END{print NR}' "${out_dir}/tmp_miss_rm_2.bim")
miss_rm=$((ct_orig - ct_miss))
echo -e "Pre-Filtering\tRemove SNPs with missingness > 5%\t()\t($miss_rm)" >> "$log_file"

# Apply MAF filter 0.000001
plink2 --bfile "${out_dir}/tmp_miss_rm_2" \
   --maf 0.000001 \
   --nonfounders \
   --make-bed \
   --out "${out_dir}/tmp_maf_rm"
check_file_exists "${out_dir}/tmp_maf_rm.bim"
ct_maf=$(awk 'END{print NR}' "${out_dir}/tmp_maf_rm.bim")
maf_rm=$((ct_miss - ct_maf))
echo -e "Pre-Filtering\tRemove SNPs with MAF < 0.000001\t()\t($maf_rm)" >> "$log_file"

# Make non-tmp version of files
cp "${out_dir}/tmp_maf_rm.bim" "${out_dir}/pre_qc_prefilt.bim"
cp "${out_dir}/tmp_maf_rm.bed" "${out_dir}/pre_qc_prefilt.bed"
cp "${out_dir}/tmp_maf_rm.fam" "${out_dir}/pre_qc_prefilt.fam"

# Add section to log file
check_file_exists "${out_dir}/pre_qc_prefilt.bim"
ct_var_prefilt=$(awk 'END{print NR}' "${out_dir}/pre_qc_prefilt.bim")
ct_samp_prefilt=$(awk 'END{print NR}' "${out_dir}/pre_qc_prefilt.fam")
echo -e "Summary\tData counts after pre-filtering\t${ct_samp_prefilt}\t${ct_var_prefilt}" >> "$log_file"

#### Sex-specific -------------------------------------------------------------

# TO NOTE: this section is selected based on provided chr list, not
# based on checking what is actually available or not in input file.
if printf "%s\n" "${chr[@]}" | grep -qE '^(X)$'; then
   echo "Optional: X chromosome steps will be run."
   x_flag=1
else
   x_flag=0
fi
if printf "%s\n" "${chr[@]}" | grep -qE '^(Y)$'; then
   echo "Optional: Y chromosome steps will be run."
   y_flag=1
else
   y_flag=0
fi

# Sex-specific functions
# Remove Y-chr SNPs with given missingness in males
filt_chrY_male_miss() {
   local bfile="$1"
   local geno="$2"
   local bfile_out="$3"

   plink2 --bfile "${bfile}" \
      --filter-males \
      --chr Y \
      --geno ${geno} \
      --nonfounders \
      --write-snplist \
      --out "${out_dir}/tmp_chrY_male_keep"
   plink2 --bfile "${bfile}" \
      --filter-males \
      --chr Y \
      --nonfounders \
      --write-snplist \
      --out "${out_dir}/tmp_chrY_male_all"

   comm -23 <(sort "${out_dir}/tmp_chrY_male_all.snplist") \
      <(sort "${out_dir}/tmp_chrY_male_keep.snplist") > \
      "${bfile_out}.txt"

   plink2 --bfile "${bfile}" \
      --exclude "${bfile_out}.txt" \
      --make-bed \
      --out "${bfile_out}"

   check_file_exists "${bfile_out}.txt"
   local miss_rm=$(awk 'END{print NR}' "${bfile_out}.txt")
   local geno_perc=$(echo "scale=0; $geno*100/1" | bc)
   echo -e "Sex-Specific\tRemove Y-chr SNPs with missingness >$geno_perc% in males\t()\t($miss_rm)" >> "$log_file"
}

# Remove X-chr SNPs with given heterozygosity in males
filt_chrX_male_het() {
   local bfile="$1"
   local perc_het="$2"
   local bfile_out="$3"

   plink2 --bfile "${bfile}" \
      --chr X \
      --split-par "$orig_build" \
      --filter-males \
      --nonfounders \
      --make-bed \
      --out "${bfile_out}_subset"

   # Pass through PLINK1.9 to generate .hh file
   plink --bfile "${bfile_out}_subset" \
      --keep-allele-order \
      --make-bed \
      --out "${bfile_out}"

   check_file_exists "${bfile_out}.hh"
   check_file_exists "${bfile_out}.fam"
   n_male=$(awk 'END{print NR}' "${bfile_out}.fam")
   Rscript ${code_dir}/scripts/filt_chrX_het.R \
      -i "${bfile_out}.hh" \
      -m "$n_male" \
      -p "$perc_het" \
      -l "$log_file"

   plink2 --bfile "${bfile}" \
      --exclude "${bfile_out}.txt" \
      --make-bed \
      --out "${bfile_out}"
}

# Remove Y-chr SNPs with genotypes at given percentage in females
filt_chrY_female_gt() {
   local bfile="$1"
   local geno="$2"
   local bfile_out="$3"

   plink2 --bfile "${bfile}" \
      --chr Y \
      --filter-females \
      --geno ${geno} \
      --nonfounders \
      --write-snplist \
      --out "${out_dir}/tmp_chrY_female_keep"

   plink2 --bfile "${bfile}" \
      --chr Y \
      --filter-females \
      --nonfounders \
      --write-snplist \
      --out "${out_dir}/tmp_chrY_female_all"

   comm -23 <(sort "${out_dir}/tmp_chrY_female_all.snplist") \
      <(sort "${out_dir}/tmp_chrY_female_keep.snplist") > \
      "${bfile_out}.txt"

   plink2 --bfile ${bfile} \
      --exclude "${bfile_out}.txt" \
      --make-bed \
      --out "${bfile_out}"

   check_file_exists "${bfile_out}.txt"
   local miss_rm=$(awk 'END{print NR}' "${bfile_out}.txt")
   local geno_perc=$(echo "scale=0; $geno*100/1" | bc)
   echo -e "Sex-Specific\tRemove Y-chr SNPs with genotypes >$geno_perc% in females\t()\t($miss_rm)" >> "$log_file"
}

# Remove Y-chr SNPs with missing >20% in males
if [ "$y_flag" -eq 1 ]; then
   echo "Optional: filtering chrY SNPs with missingness >20% in males."
   filt_chrY_male_miss \
      "${out_dir}/pre_qc_prefilt" \
      0.2 \
      "${out_dir}/tmp_chrY_male_miss_rm"
else
   cp "${out_dir}/pre_qc_prefilt.bim" "${out_dir}/tmp_chrY_male_miss_rm.bim"
   cp "${out_dir}/pre_qc_prefilt.bed" "${out_dir}/tmp_chrY_male_miss_rm.bed"
   cp "${out_dir}/pre_qc_prefilt.fam" "${out_dir}/tmp_chrY_male_miss_rm.fam"
fi

# Remove X-chr SNPs with heterozygosity >5% in males
if [ "$x_flag" -eq 1 ]; then
   echo "Optional: filtering chrX SNPs with heterozygosity >5% in males."
   filt_chrX_male_het \
      "${out_dir}/tmp_chrY_male_miss_rm" \
      5 \
      "${out_dir}/tmp_chrX_het_male_rm"
else
   cp ${out_dir}/tmp_chrY_male_miss_rm.bim ${out_dir}/tmp_chrX_het_male_rm.bim
   cp ${out_dir}/tmp_chrY_male_miss_rm.bed ${out_dir}/tmp_chrX_het_male_rm.bed
   cp ${out_dir}/tmp_chrY_male_miss_rm.fam ${out_dir}/tmp_chrX_het_male_rm.fam
fi

# Remove Y-chr SNPs with genotypes >10% in females
if [ "$y_flag" -eq 1 ]; then
   echo "Optional: filtering chrY SNPs with genotypes >10% in females."
   filt_chrY_female_gt \
      "${out_dir}/tmp_chrX_het_male_rm" \
      0.1 \
      "${out_dir}/tmp_chrY_female_miss_rm"
else
   cp "${out_dir}/tmp_chrX_het_male_rm.bim" "${out_dir}/tmp_chrY_female_miss_rm.bim"
   cp "${out_dir}/tmp_chrX_het_male_rm.bed" "${out_dir}/tmp_chrY_female_miss_rm.bed"
   cp "${out_dir}/tmp_chrX_het_male_rm.fam" "${out_dir}/tmp_chrY_female_miss_rm.fam"
fi

# If have chrX, remove mislabeled as male & female
# TO NOTE: only using chrX to check sex, could update in future
if [ "$x_flag" -eq 1 ]; then
   echo "Optional: removing samples with mislabeled sex based on chrX."
   echo "Warning: the plot sexcheck_plot.png should be used to adjust the "
   echo "default thresholds based on each dataset's distribution."

   plink2 --bfile "${out_dir}/tmp_chrY_female_miss_rm" \
      --split-par "$orig_build" \
      --check-sex min-male-xf=${sexcheck[0]} max-female-xf=${sexcheck[1]} \
      --nonfounders \
      --out "${out_dir}/tmp_sexcheck"

   check_file_exists "${out_dir}/tmp_sexcheck.sexcheck"
   Rscript ${code_dir}/scripts/sex_check.R \
      -s "${out_dir}/tmp_sexcheck.sexcheck" \
      -m "${sexcheck[0]}" \
      -f "${sexcheck[1]}" \
      -l "$log_file"

   plink2 --bfile "${out_dir}/tmp_chrY_female_miss_rm" \
      --remove "${out_dir}/tmp_sexcheck_rm.txt" \
      --make-bed \
      --out "${out_dir}/tmp_sexcheck_rm"
else
   cp "${out_dir}/tmp_chrY_female_miss_rm.bim" "${out_dir}/tmp_sexcheck_rm.bim"
   cp "${out_dir}/tmp_chrY_female_miss_rm.bed" "${out_dir}/tmp_sexcheck_rm.bed"
   cp "${out_dir}/tmp_chrY_female_miss_rm.fam" "${out_dir}/tmp_sexcheck_rm.fam"

fi

# Remove Y-chr SNPs with missing >5% in males
if [ "$y_flag" -eq 1 ]; then
   echo "Optional: filtering chrY SNPs with missingness >5% in males."
   filt_chrY_male_miss \
      "${out_dir}/tmp_sexcheck_rm" \
      0.05 \
      "${out_dir}/tmp_chrY_male_miss_rm_2"
else
   cp "${out_dir}/tmp_sexcheck_rm.bim" "${out_dir}/tmp_chrY_male_miss_rm_2.bim"
   cp "${out_dir}/tmp_sexcheck_rm.bed" "${out_dir}/tmp_chrY_male_miss_rm_2.bed"
   cp "${out_dir}/tmp_sexcheck_rm.fam" "${out_dir}/tmp_chrY_male_miss_rm_2.fam"
fi

# Remove X-chr SNPs with heterozygosity >1% in males
if [ "$x_flag" -eq 1 ]; then
   echo "Optional: filtering chrX SNPs with heterozygosity >1% in males."
   filt_chrX_male_het \
      "${out_dir}/tmp_chrY_male_miss_rm_2" \
      1 \
      "${out_dir}/tmp_chrX_het_male_rm_2"
else
   cp "${out_dir}/tmp_chrY_male_miss_rm_2.bim" "${out_dir}/tmp_chrX_het_male_rm_2.bim"
   cp "${out_dir}/tmp_chrY_male_miss_rm_2.bed" "${out_dir}/tmp_chrX_het_male_rm_2.bed"
   cp "${out_dir}/tmp_chrY_male_miss_rm_2.fam" "${out_dir}/tmp_chrX_het_male_rm_2.fam"
fi

# Remove Y-chr SNPs with genotypes >2% in females
if [ "$y_flag" -eq 1 ]; then
   echo "Optional: filtering chrY SNPs with genotypes >2% in females."
   filt_chrY_female_gt \
      "${out_dir}/tmp_chrX_het_male_rm_2" \
      0.02 \
      "${out_dir}/tmp_chrY_female_miss_rm_2"
else
   cp "${out_dir}/tmp_chrX_het_male_rm_2.bim" "${out_dir}/tmp_chrY_female_miss_rm_2.bim"
   cp "${out_dir}/tmp_chrX_het_male_rm_2.bed" "${out_dir}/tmp_chrY_female_miss_rm_2.bed"
   cp "${out_dir}/tmp_chrX_het_male_rm_2.fam" "${out_dir}/tmp_chrY_female_miss_rm_2.fam"
fi

if [[ "$x_flag" -eq 1 || "$y_flag" -eq 1 ]]; then
   # Make non-tmp version of files
   cp "${out_dir}/tmp_chrY_female_miss_rm_2.bim" "${out_dir}/pre_qc_sexcheck.bim"
   cp "${out_dir}/tmp_chrY_female_miss_rm_2.bed" "${out_dir}/pre_qc_sexcheck.bed"
   cp "${out_dir}/tmp_chrY_female_miss_rm_2.fam" "${out_dir}/pre_qc_sexcheck.fam"

   # Add section to log file
   check_file_exists "${out_dir}/pre_qc_sexcheck.bim"
   ct_var_sexcheck=$(awk 'END{print NR}' "${out_dir}/pre_qc_sexcheck.bim")
   ct_samp_sexcheck=$(awk 'END{print NR}' "${out_dir}/pre_qc_sexcheck.fam")
   echo -e "Summary\tData counts after sex-specific checks\t${ct_samp_sexcheck}\t${ct_var_sexcheck}" >> "$log_file"

   # Make tmp files for next section
   cp "${out_dir}/pre_qc_sexcheck.bim" "${out_dir}/tmp_pre_qc_sexcheck.bim"
   cp "${out_dir}/pre_qc_sexcheck.bed" "${out_dir}/tmp_pre_qc_sexcheck.bed"
   cp "${out_dir}/pre_qc_sexcheck.fam" "${out_dir}/tmp_pre_qc_sexcheck.fam"
else
   # Make tmp version of files in case this section not run
   cp "${out_dir}/pre_qc_prefilt.bim" "${out_dir}/tmp_pre_qc_sexcheck.bim"
   cp "${out_dir}/pre_qc_prefilt.bed" "${out_dir}/tmp_pre_qc_sexcheck.bed"
   cp "${out_dir}/pre_qc_prefilt.fam" "${out_dir}/tmp_pre_qc_sexcheck.fam"
fi

# SNP Deduplication -----------------------------------------------------------

# Deduplicate and set all variant IDs to chr:pos:ref:alt
# First, if duplicate IDs have different missingness, remove the SNP with
# more missingness.
plink2 --bfile "${out_dir}/tmp_pre_qc_sexcheck" \
   --missing --freq \
   --make-pgen \
   --out "${out_dir}/tmp_dedup"
Rscript ${code_dir}/scripts/dedup_miss.R \
   -v "${out_dir}/tmp_dedup.vmiss" \
   -p "${out_dir}/tmp_dedup.pvar" \
   -l "$log_file"
plink2 --bfile "${out_dir}/tmp_pre_qc_sexcheck" \
   --make-bed \
   --exclude "${out_dir}/tmp_dedup_rm.txt" \
   --out "${out_dir}/tmp_dedup_rm"

# Then, update all IDs to chr:pos:ref:alt and remove the remaining duplicates
# by arbitrarily keeping first
plink2 --bfile "${out_dir}/tmp_dedup_rm" \
   --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1000 \
   --rm-dup force-first \
   --make-bed \
   --out "${out_dir}/pre_qc_chrpos_dedup"

check_file_exists "${out_dir}/tmp_dedup_rm.bim"
check_file_exists "${out_dir}/pre_qc_chrpos_dedup.bim"
ct_dedup_miss=$(awk 'END{print NR}' "${out_dir}/tmp_dedup_rm.bim")
ct_dedup_arb=$(awk 'END{print NR}' "${out_dir}/pre_qc_chrpos_dedup.bim")
dedup_arb_rm=$((ct_dedup_miss - ct_dedup_arb))
echo -e "SNP-Dedup\tRemove duplicate SNPs arbitrarily\t()\t(${dedup_arb_rm})" >> "$log_file"

# Add summary section to log file
check_file_exists "${out_dir}/pre_qc_chrpos_dedup.bim"
ct_var_dedup=$(awk 'END{print NR}' "${out_dir}/pre_qc_chrpos_dedup.bim")
ct_samp_dedup=$(awk 'END{print NR}' "${out_dir}/pre_qc_chrpos_dedup.fam")
echo -e "Summary\tData counts after SNP deduplication\t${ct_samp_dedup}\t${ct_var_dedup}" >> "$log_file"

# HWE -------------------------------------------------------------------------

orig_build_num=$(echo "$orig_build" | grep -o '[0-9]\+')
to_build_num=$(echo "$to_build" | grep -o '[0-9]\+')

# Perform HWE differentially by chr6 MHC
for c in "${chr[@]}"; do
   if [ "$c" == "6" ]; then
      echo "Processing chr6 with HWE 1e-20."
      # Get chr6 MHC region from Paul Norman's coordinates
      start=$(head -n 1 "${code_dir}/refs/mhc_extended_hg${orig_build_num}.bed" |
         awk -F':' '{gsub(/-.*/, "", $2); print $2}')
      stop=$(tail -n 1 "${code_dir}/refs/mhc_extended_hg${orig_build_num}.bed" |
         awk -F':' '{gsub(/-.*/, "", $2); print $2}')

      plink2 --bfile "${out_dir}/pre_qc_chrpos_dedup" \
         --chr 6 --from-bp $start --to-bp $stop \
         --make-pgen --out "${out_dir}/tmp_mhc"

      # Get all chr6 SNPs not in region
      awk '/^#/ {next} {print $3}' "${out_dir}/tmp_mhc.pvar" > \
         "${out_dir}/chr6_mhc_var_id_list.txt"

      plink2 --bfile "${out_dir}/pre_qc_chrpos_dedup" \
         --chr 6 --exclude "${out_dir}/chr6_mhc_var_id_list.txt" \
         --make-pgen --out "${out_dir}/tmp_non_mhc"

      # Set var to track who HWE calculated in
      hwe_subset=""
      if [ -n "$hwe_filt_list" ]; then
         # Apply HWE filter only in given ID list
         hwe_subset="controls"
         echo "Optional: applying HWE filter only in given ID list."
         plink2 --pfile "${out_dir}/tmp_mhc" \
            --hwe 1e-20 --nonfounders --keep "$hwe_filt_list" \
            --write-snplist --out "${out_dir}/tmp_mhc_hwe_list"

         plink2 --pfile "${out_dir}/tmp_non_mhc" \
            --hwe 1e-6 --nonfounders --keep "$hwe_filt_list" \
            --write-snplist --out "${out_dir}/tmp_non_mhc_hwe_list"
      else
         hwe_subset="all"
         echo "Optional: applying HWE filter in all."
         plink2 --pfile "${out_dir}/tmp_mhc" \
            --hwe 1e-20 --nonfounders \
            --write-snplist --out "${out_dir}/tmp_mhc_hwe_list"

         plink2 --pfile "${out_dir}/tmp_non_mhc" \
            --hwe 1e-6 --nonfounders \
            --write-snplist --out "${out_dir}/tmp_non_mhc_hwe_list"
      fi

      plink2 --pfile "${out_dir}/tmp_mhc" \
         --extract "${out_dir}/tmp_mhc_hwe_list.snplist" \
         --make-bed --out "${out_dir}/tmp_mhc_hwe"

      plink2 --pfile "${out_dir}/tmp_non_mhc" \
         --extract "${out_dir}/tmp_non_mhc_hwe_list.snplist" \
         --make-bed --out "${out_dir}/tmp_non_mhc_hwe"
   
      # Merge chr6 back together
      plink --bfile "${out_dir}/tmp_mhc_hwe" \
         --bmerge "${out_dir}/tmp_non_mhc_hwe" \
         --keep-allele-order --allow-no-sex \
         --make-bed --out "${out_dir}/tmp_chr6"

      check_file_exists "${out_dir}/tmp_mhc.pvar"
      check_file_exists "${out_dir}/tmp_mhc_hwe.bim"
      ct_mhc=$(awk 'END{print NR}' "${out_dir}/tmp_mhc.pvar")
      ct_mhc_hwe=$(awk 'END{print NR}' "${out_dir}/tmp_mhc_hwe.bim")
      ct_mhc_hwe_rm=$((ct_mhc -1 - ct_mhc_hwe))  # -1 for .pvar header

      check_file_exists "${out_dir}/tmp_non_mhc.pvar"
      check_file_exists "${out_dir}/tmp_non_mhc_hwe.bim"
      ct_non_mhc=$(awk 'END{print NR}' "${out_dir}/tmp_non_mhc.pvar")
      ct_non_mhc_hwe=$(awk 'END{print NR}' "${out_dir}/tmp_non_mhc_hwe.bim")
      ct_non_mhc_hwe_rm=$((ct_non_mhc -1 - ct_non_mhc_hwe))  # -1 for .pvar header

      echo -e "HWE\tRemove chr6 MHC SNPs failing HWE (1e-20) in ${hwe_subset}\t()\t(${ct_mhc_hwe_rm})" >> "$log_file"
      echo -e "HWE\tRemove chr6 non-MHC SNPs failing HWE (1e-6) in ${hwe_subset}\t()\t(${ct_non_mhc_hwe_rm})" >> "$log_file"

   else
      echo "Processing chr$c with HWE 1e-6."
      plink2 --bfile "${out_dir}/pre_qc_chrpos_dedup" \
         --chr "$c" \
         --make-pgen --out "${out_dir}/tmp_chr${c}_no_hwe"

      if [ -n "$hwe_filt_list" ]; then
         # Apply HWE filter only in given ID list
         echo "Optional: applying HWE filter only in given ID list."
         plink2 --pfile "${out_dir}/tmp_chr${c}_no_hwe" \
            --hwe 1e-6 --nonfounders --keep "$hwe_filt_list" \
            --write-snplist --out "${out_dir}/tmp_chr${c}_hwe_list"
      else
         echo "Optional: applying HWE filter in all."
         plink2 --pfile "${out_dir}/tmp_chr${c}_no_hwe" \
            --hwe 1e-6 --nonfounders \
            --write-snplist --out "${out_dir}/tmp_chr${c}_hwe_list"
      fi

      plink2 --pfile "${out_dir}/tmp_chr${c}_no_hwe" \
         --extract "${out_dir}/tmp_chr${c}_hwe_list.snplist" \
         --make-bed --out "${out_dir}/tmp_chr${c}"
   fi
done

# If multiple chrs being processed, prep for merge of all chr,
# which is needed to next pipeline step
if [ "${#chr[@]}" -gt 1 ]; then
   echo "Processing more than one chr."
   base_chr=${chr[0]}
   merge_list="${out_dir}/tmp_merge_list.txt"

   for c in "${chr[@]}"; do
      if [ "$c" != "$base_chr" ]; then
         echo "${out_dir}/tmp_chr${c}" >> "$merge_list"
      fi
   done

   # Merge all chr currently preparing
   plink --bfile "${out_dir}/tmp_chr${base_chr}" \
      --merge-list "$merge_list" \
      --keep-allele-order --allow-no-sex --output-chr MT \
      --make-bed --out "${out_dir}/pre_qc_hwe"
   
else
   echo "Processing only one chr."
   plink --bfile "${out_dir}/tmp_chr${chr[0]}" \
      --keep-allele-order --output-chr MT \
      --make-bed --out "${out_dir}/pre_qc_hwe"
fi

# Add HWE for all chr preparing to log file
check_file_exists "${out_dir}/pre_qc_hwe.bim"
ct_var_hwe=$(awk 'END{print NR}' "${out_dir}/pre_qc_hwe.bim")
ct_samp_hwe=$(awk 'END{print NR}' "${out_dir}/pre_qc_hwe.fam")
ct_hwe_rm=$((ct_var_dedup - ct_var_hwe))
echo -e "HWE\tRemove SNPs (all chr) failing HWE (1e-20 MHC, 1e-6 otherwise) in ${hwe_subset}\t()\t(${ct_hwe_rm})" >> "$log_file"
echo -e "Summary\tData counts after HWE\t${ct_samp_hwe}\t${ct_var_hwe}" >> "$log_file"

# Liftover --------------------------------------------------------------------

# Check if liftover needed
if [ "$orig_build_num" != "$to_build_num" ]; then
   echo "Optional: Lifting over from ${orig_build} to ${to_build}."

   # Create bed file to crossover from hg19 to hg38 
   cat "${out_dir}/pre_qc_hwe.bim" | cut -f1 | sed 's/^/chr/' > "${out_dir}/tmp_c1.txt"
   cat "${out_dir}/pre_qc_hwe.bim" | cut -f4 > "${out_dir}/tmp_c2.txt"
   cat "${out_dir}/pre_qc_hwe.bim" | cut -f4 > "${out_dir}/tmp_c3.txt"
   cat "${out_dir}/pre_qc_hwe.bim" | cut -f2 > "${out_dir}/tmp_c4.txt"
   paste "${out_dir}/tmp_c1.txt" \
         "${out_dir}/tmp_c2.txt" \
         "${out_dir}/tmp_c3.txt" \
         "${out_dir}/tmp_c4.txt" \
         >  "${out_dir}/tmp_in.bed"

   CrossMap bed "${code_dir}/refs/hg${orig_build_num}ToHg${to_build_num}.over.chain" \
      "${out_dir}/tmp_in.bed"  \
      "${out_dir}/tmp_out.bed"

   # Extract only those SNPs that were successfully cross-overed
   cut -f4 "${out_dir}/tmp_out.bed" > "${out_dir}/tmp_snp_keep.txt"
   plink --bfile "${out_dir}/pre_qc_hwe" \
      --extract "${out_dir}/tmp_snp_keep.txt" \
      --keep-allele-order \
      --make-bed --out "${out_dir}/pre_qc_hg${to_build_num}"

   # Update bim file positions
   Rscript --vanilla ${code_dir}/scripts/update_pos.R \
      -c "${out_dir}/tmp_out.bed" -s "${out_dir}/pre_qc_hg${to_build_num}.bim"

   # Add liftover summary to log file
   check_file_exists "${out_dir}/pre_qc_hg${to_build_num}.bim"
   ct_var_liftover=$(awk 'END{print NR}' "${out_dir}/pre_qc_hg${to_build_num}.bim")
   ct_liftover_rm=$((ct_var_hwe - ct_var_liftover))
   echo -e "Imputation-Prep\tRemove SNPs failing liftover to hg${to_build_num}\t()\t(${ct_liftover_rm})" >> "$log_file"

   # Make tmp files
   cp "${out_dir}/pre_qc_hg${to_build_num}.bim" "${out_dir}/tmp_liftover.bim"
   cp "${out_dir}/pre_qc_hg${to_build_num}.bed" "${out_dir}/tmp_liftover.bed"
   cp "${out_dir}/pre_qc_hg${to_build_num}.fam" "${out_dir}/tmp_liftover.fam"
else 
   echo "Optional: Not lifting over."
   cp "${out_dir}/pre_qc_hwe.bim" "${out_dir}/tmp_liftover.bim"
   cp "${out_dir}/pre_qc_hwe.bed" "${out_dir}/tmp_liftover.bed"
   cp "${out_dir}/pre_qc_hwe.fam" "${out_dir}/tmp_liftover.fam"
fi

# Imputation prep -------------------------------------------------------------

# Remove strand ambiguous SNPs, set IDs to chr:pos:ref:alt
check_file_exists "${out_dir}/tmp_liftover.bim"
Rscript --vanilla ${code_dir}/scripts/get_strand_amb_SNPs.R \
   -b "${out_dir}/tmp_liftover.bim"
plink2 --bfile "${out_dir}/tmp_liftover" \
   --exclude "${out_dir}/tmp_strand_remove_snps.txt" \
   --make-bed --out "${out_dir}/tmp_pre_qc"
plink2 --bfile "${out_dir}/tmp_pre_qc" \
   --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 1000 \
   --make-bed --out "${out_dir}/pre_qc"

# Write out VCF files split by chr for imputation.
for c in "${chr[@]}"; do
   plink --bfile "${out_dir}/pre_qc" \
      --chr $c --keep-allele-order \
      --recode vcf --out "${out_dir}/tmp_chr${c}_exp"
   if [ "$to_build_num" == "38" ]; then
      vcf-sort "${out_dir}/tmp_chr${c}_exp.vcf" | \
         sed -E 's/^([0-9XYM]+)/chr\1/' | \
         bgzip -c > "${out_dir}/chr${c}_pre_qc.vcf.gz"
   else
      vcf-sort "${out_dir}/tmp_chr${c}_exp.vcf" | \
         sed -E 's/^chr([0-9XYM]+)/\1/' | \
         bgzip -c > "${out_dir}/chr${c}_pre_qc.vcf.gz"
   fi
done

# Add ambiguous SNPs and final data counts
check_file_exists "${out_dir}/pre_qc.bim"
ct_var_final=$(awk 'END{print NR}' "${out_dir}/pre_qc.bim")
ct_samp_final=$(awk 'END{print NR}' "${out_dir}/pre_qc.fam")
ct_before_ambig=$(awk 'END{print NR}' "${out_dir}/tmp_liftover.bim")
ct_ambig_rm=$((ct_before_ambig - ct_var_final))
echo -e "Imputation-Prep\tRemove ambiguous A/T C/G SNPs\t()\t(${ct_ambig_rm})" >> "$log_file"
echo -e "Summary\tFinal data counts\t${ct_samp_final}\t${ct_var_final}" >> "$log_file"

# Cleanup
rm ${out_dir}/tmp_*
