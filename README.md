## TODOs

* update below!
* add instructions for getting TOPMed API key and/or manual submit?
* add note that need to lift input data over to match reference panel
    build otherwise won't be able to use fix strands code.
* add checks that config values are only the allowable ones?
* add more documentation to python functions
* think could remove a config arg so that to_build is auto set based on chosen
    imputation
* not sure how to better handle relative paths within the pipeline directory?
* add documentation about geno in all individuals followed by hwe in controls -
    may lead to warning since subsetting increases variant missingness, but
    think should be okay to ignore: https://groups.google.com/g/plink2-users/c/xjmPNfN0Swc/m/79_CAtMRBwAJ
* better handle using imputationbot v1 for TOPmed and v2 for Michigan
    * https://github.com/UW-GAC/primed-imputation/blob/main/register_token.sh
    
    
## TODOs - after 10/29/25
* add R argparse and aria2 to mamba env
* discuss sex chromosome options, note that SEX column in PLINK file must be filled out
* add PLINK2 to conda recipe - but needs to be OSX64 or linux?
* TODO in future - convert this to Nextflow?
* awknowledge that better than duplicate removal would be to consider fixin mismapped
    SNPs using reference panel... potential future step
* add note that create initial input should be run twice so that plot with sex checks
    can be examined and thresholds can be adjusted
* get job-id from log files
* update create_initial_input log file to count all non-chr6 HWE (vs all chr HWE)
* add note the pre-filtering, sex-checks done for all chr present in given file, but
    SNP dedup, HWE and imputation prep only done for requested chr to prepare
* add note that chrY can be prepared but not imputed
* add concat_convert_to_plink as main option
* split into two snakefiles - one pre- and one post-QC (would help with post QC paths
    being different)
* need to remove imputed files that don't want to keep
* need to confirm download works for Michigan imputation


## **imputation_snakemake**

Pipeline for imputing autosomal genetic array data with hg19 or hg38 coordinates
to the TOPMed r3 reference panel, the Michigian 1000 genomes phase 3 version 5
reference panel, or the Michigan HLA four-digit multi-ethnic v1 or v2 panels.

Adapted from Michelle Daya's topmed_imputation pipeline.

## Setup

### Apptainer

To build apptainer on an x86_64 operating system, from the top level of this
repository run:

``` 
apptainer build envs/topmed_imputation.sif envs/topmed_imputation.def

```

### Conda Environment

TODO: add instructions to use the recipe file at envs/env_imputation_snakemake.yml

## Input Files

* Genetic array data in PLINK1.9 file format.

## How to Run

### Step 1

Copy "config.yml" to the project repo want to run imputation for, and update paths and
settings specific to the data:

```
### TODO: Update pipeline settings
# Default filters of MAF 0 and Rsq 0.3 are applied after imputation.
chr: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]  # list of chromosome numbers with format [1,10,22]
orig_build: "hg19"  #  either "hg19" or "hg38"
to_build: "hg38"  # either "hg19" or "hg38"
imp: "topmed"  # should be 'topmed', 'mich_hla_v1', 'mich_hla_v2', 'mich_1kg_p3_v5', 'mich_hrc'
imp_rsq_filt: "0"  # filter to apply on imputation server (0, 0.001, 0.1, 0.2, or 0.3)
imp_name: "job-name"  # job name for imputation server
imp_job_id: ""  # should be "job-####-##-###"
zip_pw: ""  # add from email when imputed job finishes
opt: "gt"  # option for "gt" only or "all" dosage information (much slower)
use_cont: false  # true to use container, false to use local machine

### TODO: Update host paths outside container
# PLINK1.9 file prefix with path - must be single file, but can contain one or more chr
plink_prefix: "/path/to/plink"
# ID list to calculate HWE filtering in (default is controls)
id_list_hwe: "/path/to/crtls_id_list.txt"
# Top-level directory for all output
out_dir: "/path/to/output"
# Repo with snakemake pipeline
repo: "/Users/slacksa/repos/imputation_snakemake"

```

Note that some variables (imp_job_id & zip_pw) will need to be added to the config file as
the pipeline is run.

Note that imp_rsq_filt refers to a filter based on imputation quality Rsq that would be
applied by the imputation server. To preserve information about all variants the server
attempted to impute, it is best to set this to 0. If necessary, due to large file size,
a filter can be applied on the server by setting this option. However, although the server
generally returns TYPED input variants, it does almost attempt to impute them and applying
the server-side Rsq filter appears to removed TYPED variants that imputed with an Rsq
less than the filter as of 10/29/25.

Note that for TOPMed imputation, the 'population' option is set to 'all', which means allele
frequencies will be compared between the input data and the TOPMed panel to generate a QC
report. For Michigan imputation, the 'population' option is set to 'off', since (at least
for 1000G whole genome imputation) there is no good option for allele frequency comparison
for mixed populations. Either way, imputation results should not be affected.

*TODO: add note about activating conda environment if set use_cont to false.*

### Step 1

Use the "run_snakemake.sh" script to run the snakefile through initial data preparation and
submission to the imputation server for QC only. This will create and upload the initial
input files to the imputation server and panel selected in the config file.

```
bash run_snakemake.sh \
    -s submit_initital_input \
    -f /path/to/project/config.yml

```

Note that you can add "-d" to run a snakemake --dry-run, use "-c" to set the number of
cores used (the default is 6), or use "-h" to show the run script help message.

This will create and upload the initial input files to the imputation server and panel
selected in the config file.

The pre-imputation QC includes liftover (if necessary), removal of strand ambiguous SNPs,
updating all variant IDs to chr:pos:ref:alt, filtering by PLINK2 --maf (1e-6), --geno (0.05),
and --hwe (1e-20 for chr6 MHC region, 1e-6 for all other regions/chromosomes). A summary of
each step is included in the log file.

## Step 2

*TODO: is the auto-download part working?*

Once the QC automatically submitted has run, you will receive an email. The pipeline here
will download the results for you, modify the input based on the server QC, and resubmit
to the server for imputation. However, you must log onto the imputation server, get the
job ID, and add it to the project-specific config.yml first:

```
### TODO: Update pipeline settings
# ...
imp_job_id: "job-1234-56-789"  # should be "job-####-##-###"
zip_pw: ""  # add from email when imputed job finishes

```

Only add the imp_job_id for now. You will add the zip_pw later.

Use the "run_snakemake.sh" to run the pipeline through the "submit_fix_strands" step:

```
bash run_snakemake.sh \
    -s submit_fix_strands \
    -f /path/to/project/config.yml

```

This step flips strands of variants identified as such in the snps-exlcuded.txt file - this
will produce the final post QC VCF files for imputation.

*TO NOTE: if using Michigan imputation server, may need to move files out of*
*subfolders manually.*

## Step 3

Once the imputation has run, you will receive an email. The pipeline here
will download the results for you, filter them, and create PLINK files. However,
you must log onto the imputation server, get the updated job ID, and add it and 
the zip password that was emailed to you to the project-specific config.yml first:

```
### TODO: Update pipeline settings
# ...
imp_job_id: "job-2345-67-891"  # should be "job-####-##-###"
zip_pw: "abc123zyx987"  # add from email when imputed job finishes

```

Use the "run_snakemake.sh" to run the pipeline through the "filter_info_and_vcf_files" step:

```
bash run_snakemake.sh \
    -s filter_info_and_vcf_files \
    -f /path/to/project/config.yml

```
