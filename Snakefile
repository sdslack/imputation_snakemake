# Setup -----------------------------------------------------------------------------------------

### Get variables from config file
# Required pipeline settings
chr: List[str] = config["chr"]
orig_build: str = config["orig_build"]
to_build: str = config["to_build"]
imp: str = config["imp"]
imp_name: str = config["imp_name"]
zip_pw: str = config["zip_pw"]
imp_job_id: str = config["imp_job_id"]

# Required pipeline settings with defaults
imp_rsq_filt: str = config.get("imp_rsq_filt", "0")
opt: str = config.get("opt", "gt")
use_cont: bool = config.get("use_cont", True)
sexcheck: List[str] = config.get("sexcheck", ["0.8", "0.2"])

# Required host paths outside container
plink_prefix: str = config["plink_prefix"]
plink_prefix_name = Path(plink_prefix).name

# Optional host paths outside container
id_list: str = config.get("id_list", None)
id_list_hwe: str = config.get("id_list_hwe", None)

if use_cont:
    # Container paths
    plink_dir: str = config["plink_dir_cont"]
    out_dir: str = config["out_dir_cont"]
else:
    plink_dir = Path(plink_prefix).parent
    out_dir: str = config["out_dir"]

# Modify arrays for bash script
chr_str = " ".join(map(str, chr))
chr_noY = [c for c in chr if c != "Y"]
chr_noY_str = " ".join(map(str, chr_noY))
sexcheck_str = " ".join(map(str, sexcheck))

### Other prep
# Set default values currently not controlled by config file
maf = "0"
rsq = "0.3"

# Get dir of pipeline
code_dir = Path(workflow.snakefile).resolve().parent

# Rules -----------------------------------------------------------------------------------------

rule all:
    input:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pvar",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.psam",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pgen"

rule create_initial_input:
    input:
        f"{plink_dir}/{plink_prefix_name}.bed",
        f"{plink_dir}/{plink_prefix_name}.bim",
        f"{plink_dir}/{plink_prefix_name}.fam"
    output:
        [f"{out_dir}/pre_qc/chr{c}_pre_qc.vcf.gz" for c in chr]
    log:
        f"{out_dir}/pre_qc/create_initial_input.log"
    params:
        script=Path(code_dir, "scripts/create_initial_input.sh"),
        summary=f"{out_dir}/pre_qc/summary_pre_qc.txt"
    run:
        cmd = (
            f"bash {params.script} "
            f"-p {plink_dir}/{plink_prefix_name} "
            f"-o {out_dir}/pre_qc "
            f"-c {code_dir} "
            f"-n \"{chr_str}\" "
            f"-b {orig_build} "
            f"-t {to_build} "
            f"-s \"{sexcheck_str}\" "
            f"-l {params.summary} "
        )
        # If id_list given for filtering
        if id_list:
            id_list_name=Path(id_list).name
            if use_cont:
                id_list_dir: str = config["id_list_dir_cont"]
            else:
                id_list_dir = Path(id_list).parent
            cmd += f" -k {id_list_dir}/{id_list_name}"
        # If id_list_hwe given for HWE calculation
        if id_list_hwe:
            id_list_hwe_name=Path(id_list_hwe).name
            if use_cont:
                id_list_hwe_dir: str = config["id_list_dir_cont"]
            else:
                id_list_hwe_dir = Path(id_list_hwe).parent
            cmd += f" -h {id_list_hwe_dir}/{id_list_hwe_name}"
        # Redirect logs
        cmd += f" > {log} 2>&1"
        # Run
        shell(cmd)

# Need "" around chr_str to keep all as single input
rule submit_initial_input:
    input:
        [f"{out_dir}/pre_qc/chr{c}_pre_qc.vcf.gz" for c in chr]
    output:
        log_final=f"{out_dir}/pre_qc/submit_initial_input.log"
    params:
        script=Path(code_dir, "scripts/submit.py"),
        log_tmp=f"{out_dir}/pre_qc/tmp_submit_initial_input.log"
    shell:
        """
        python {params.script} \
            --dir {out_dir}/pre_qc \
            --chr "{chr_noY_str}" \
            --imp {imp} \
            --build {to_build} \
            --mode "qconly" \
            --imp-name {imp_name} \
            > {params.log_tmp} 2>&1

        mv {params.log_tmp} {output.log_final}
        """

rule download_qc:
    output:
        log_final=f"{out_dir}/pre_qc/download_qc.log",
        snps_excl=f"{out_dir}/pre_qc/snps-excluded.txt"
    params:
        script=Path(code_dir, "scripts/download_results.sh"),
        log_tmp=f"{out_dir}/pre_qc/tmp_download_qc.log"
    shell:
        """
        bash {params.script} \
            -i {imp} \
            -c {code_dir} \
            -o {out_dir}/pre_qc \
            -j {imp_job_id} \
            > {params.log_tmp} 2>&1

        mv {params.log_tmp} {output.log_final}
        """

rule fix_strands:
    input:
        f"{out_dir}/pre_qc/snps-excluded.txt",
        f"{out_dir}/pre_qc/download_qc.log"
    output:
        [f"{out_dir}/post_qc/chr{c}_post_qc.vcf.gz" for c in chr_noY]
    log:
        f"{out_dir}/post_qc/fix_strands.log"
    params:
        script=Path(code_dir, "scripts/fix_strands.sh")
    shell:
        """
        bash {params.script} \
            -o {out_dir} \
            -c {code_dir} \
            -n "{chr_noY_str}" \
            -t {to_build} \
            -i {imp} \
            > {log} 2>&1
        """

rule submit_fix_strands:
    input:
        [f"{out_dir}/post_qc/chr{c}_post_qc.vcf.gz" for c in chr_noY]
    output:
        log_final=f"{out_dir}/post_qc/submit_fix_strands.log"
    params:
        script=Path(code_dir, "scripts/submit.py"),
        log_tmp=f"{out_dir}/post_qc/tmp_submit_fix_strands.log"
    shell:
        """
        python {params.script} \
            --dir {out_dir}/post_qc \
            --chr "{chr_noY_str}" \
            --imp {imp} \
            --build {to_build} \
            --mode imputation \
            --rsq-filt {imp_rsq_filt} \
            --imp-name {imp_name} \
            > {params.log_tmp} 2>&1

        mv {params.log_tmp} {output.log_final}
        """

rule download_results:
    threads: 16
    output:
        [f"{out_dir}/imputed/chr_{c}.zip" for c in chr_noY]
    log:
        f"{out_dir}/imputed/download_results.log"
    params:
        script=Path(code_dir, "scripts/download_results.sh")
    shell:
        """
        bash {params.script} \
            -i {imp} \
            -c {code_dir} \
            -o {out_dir}/imputed \
            -j {imp_job_id} \
            > {log} 2>&1
        """

# Different from the other rules, this script in this rule runs once for each chr
rule unzip_results:
    input:
        [f"{out_dir}/imputed/chr{c}.dose.vcf.gz" for c in chr_noY]

rule unzip_results_helper:
    input:
        f"{out_dir}/imputed/chr_{{chr}}.zip"
    output:
        f"{out_dir}/imputed/chr{{chr}}.dose.vcf.gz"
    log:
        f"{out_dir}/imputed/chr{{chr}}_unzip_results.log"
    params:
        script=Path(code_dir, "scripts/unzip_results.sh")
    shell:
        """
        bash {params.script} \
            -d {out_dir}/imputed \
            -p "{zip_pw}" \
            -c {wildcards.chr} \
            > {log} 2>&1
        """

# Different from the other rules, this script in this rule runs once for each chr
rule filter_info_and_vcf_files:
    input:
        [f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{c}_clean.vcf.gz" for c in chr_noY]

rule filter_info_and_vcf_files_helper:
    input:
        f"{out_dir}/imputed/chr{{chr}}.dose.vcf.gz"
    output:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{{chr}}_clean.vcf.gz"
    log:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{{chr}}_filter_info_and_vcf_files.log"
    params:
        script=f'{code_dir}/scripts/filter_info_and_vcf_files{"_bcftools" if opt == "all" else ""}.sh'
    shell:
        """
        bash {params.script} \
            -n {wildcards.chr} \
            -r {rsq} \
            -m {maf} \
            -p {out_dir}/pre_qc/pre_qc \
            -d {out_dir} \
            -o {opt} \
           > {log} 2>&1
        """

rule concat_convert_to_plink:
    input:
        [f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{c}_clean.vcf.gz" for c in chr_noY]
    output:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pvar",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.psam",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pgen"
    log:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/concat_convert_to_plink.log"
    params:
        script=Path(code_dir, "scripts/concat_convert_to_plink.sh")
    shell:
        """
        bash {params.script} \
            -p {out_dir}/pre_qc/pre_qc \
            -d {out_dir}/imputed_clean_maf{maf}_rsq{rsq} \
            > {log} 2>&1
        """
