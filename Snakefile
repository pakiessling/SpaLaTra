import glob
import os

configfile: "config/config.yaml"

_input = config["input"]
if os.path.isfile(_input) and _input.endswith(".h5ad"):
    INPUT_FILES = [_input]
    SAMPLE_IDS = [os.path.splitext(os.path.basename(_input))[0]]
    INPUT_DIR = os.path.dirname(os.path.abspath(_input))
else:
    INPUT_FILES = glob.glob(os.path.join(_input, '*.h5ad'))
    SAMPLE_IDS = [os.path.splitext(os.path.basename(f))[0] for f in INPUT_FILES]
    INPUT_DIR = _input

def get_sample_input(wildcards):
    if os.path.isfile(config["input"]) and config["input"].endswith(".h5ad"):
        return config["input"]
    return os.path.join(INPUT_DIR, wildcards.sample + ".h5ad")

ALL_METHODS = ["tacco", "singler", "rctd", "phispace", "insitutype"]
METHODS = config.get("methods", ALL_METHODS)

for m in METHODS:
    if m not in ALL_METHODS:
        raise ValueError(f"Unknown method '{m}'. Valid options: {ALL_METHODS}")
if len(METHODS) < 2:
    raise ValueError("At least 2 methods must be enabled for a meaningful consensus.")

for m in METHODS:
    os.makedirs(os.path.join(config["output"], m), exist_ok=True)

LOG_DIR = config.get("log_dir", os.path.join(config["output"], "logs"))
os.makedirs(LOG_DIR, exist_ok=True)

wildcard_constraints:
    sample  = r"[^/]+",
    section = r"[^/]+"

def method_outputs(methods, sample_ids, output_dir):
    return [
        expand(os.path.join(output_dir, m, "{sample}_" + m + ".csv"), sample=sample_ids)
        for m in methods
    ]

localrules: all

rule all:
    input:
        method_outputs(METHODS, SAMPLE_IDS, config["output"]),
        os.path.join(config["output"], "consensus.csv"),
        os.path.join(config["output"], "report.html")
    shell:
        """
        echo Annotation complete!
        """

rule singler:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "singler", "{sample}_singler.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(LOG_DIR,"singler_{sample}.log")
    resources:
        mem_mb=120000,
        cpus_per_task=20
    shell:
        "Rscript scripts/run_singler.R --input {input} --ref {config[ref]} --ref_column {config[ref_column]} --output {output} > {log} 2>&1"

rule tacco:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "tacco", "{sample}_tacco.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(LOG_DIR,"tacco_{sample}.log")
    resources:
        mem_mb=100_000,
        cpus_per_task=4
    params:
        sample_col_arg = lambda wildcards: (
            f"--sample_column {config['sample_column']}"
            if config.get('sample_column') else ""
        )
    shell:
        "python scripts/run_tacco.py --input {input} --ref {config[ref]} --ct_column {config[ref_column]} --output {output} {params.sample_col_arg} > {log} 2>&1"

rule phispace:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "phispace", "{sample}_phispace.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(LOG_DIR,"phispace_{sample}.log")
    resources:
        mem_mb=150000,
        cpus_per_task=4
    shell:
        "Rscript scripts/run_phi_space.R --input {input} --ref {config[ref]} --ref_column {config[ref_column]} --output {output} > {log} 2>&1"

if config.get("sample_column"):

    def _load_section_manifest(path):
        mapping = {}
        with open(path) as fh:
            for line in fh:
                sane, orig = line.rstrip("\n").split("\t", 1)
                mapping[sane] = orig
        return mapping

    checkpoint rctd_list_sections:
        input:  get_sample_input
        output: os.path.join(config["output"], "rctd", "sections", "{sample}_sections.txt")
        log:    os.path.join(LOG_DIR, "rctd_list_sections_{sample}.log")
        conda:  "environment.yml"
        resources: mem_mb=8_000, cpus_per_task=1
        shell:
            "python scripts/list_rctd_sections.py "
            "--input {input} "
            "--sample_column {config[sample_column]} "
            "--output {output} > {log} 2>&1"

    rule rctd_one_section:
        input:
            query    = get_sample_input,
            manifest = lambda wc: checkpoints.rctd_list_sections.get(
                           sample=wc.sample).output[0]
        output:
            os.path.join(config["output"], "rctd", "sections", "{sample}", "{section}_rctd.csv")
        log:
            os.path.join(LOG_DIR, "rctd_{sample}_{section}.log")
        conda:  "environment.yml"
        threads: 5
        resources: mem_mb=50_000, cpus_per_task=5
        params:
            orig_section = lambda wc: _load_section_manifest(
                checkpoints.rctd_list_sections.get(sample=wc.sample).output[0]
            )[wc.section],
            ref_column = config.get("ref_column", "cell_subtype")
        shell:
            "Rscript scripts/run_rctd.R "
            "--input {input.query} --ref {config[ref]} "
            "--ref_column {params.ref_column} --output {output} "
            "--max_cores {threads} "
            "--sample_column {config[sample_column]} "
            "--section '{params.orig_section}' "
            "> {log} 2>&1"

    def _rctd_section_csvs(wildcards):
        manifest = checkpoints.rctd_list_sections.get(sample=wildcards.sample).output[0]
        sane_names = list(_load_section_manifest(manifest).keys())
        return expand(
            os.path.join(config["output"], "rctd", "sections", "{sample}", "{section}_rctd.csv"),
            sample=wildcards.sample, section=sane_names
        )

    rule rctd_merge_sections:
        input:  _rctd_section_csvs
        output: os.path.join(config["output"], "rctd", "{sample}_rctd.csv")
        log:    os.path.join(LOG_DIR, "rctd_merge_{sample}.log")
        conda:  "environment.yml"
        resources: mem_mb=16_000, cpus_per_task=1
        run:
            import pandas as pd
            dfs = [pd.read_csv(f, index_col=0) for f in input]
            pd.concat(dfs).to_csv(output[0])

else:
    rule rctd:
        input:
            get_sample_input
        output:
            os.path.join(config["output"], "rctd", "{sample}_rctd.csv")
        conda:
            "environment.yml"
        log:
            os.path.join(LOG_DIR,"rctd_{sample}.log")
        threads: 5
        resources:
            mem_mb=150000,
            cpus_per_task=5
        shell:
            "Rscript scripts/run_rctd.R --input {input} --ref {config[ref]} "
            "--ref_column {config[ref_column]} --output {output} "
            "--max_cores {threads} > {log} 2>&1"

rule insitutype:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "insitutype", "{sample}_insitutype.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(LOG_DIR,"insitutype_{sample}.log")
    resources:
        mem_mb=50000,
        cpus_per_task=4
    shell:
        "Rscript scripts/run_insitutype.R --input {input} --ref {config[ref]} --ref_column {config[ref_column]} --output {output} > {log} 2>&1"

rule consensus:
    input:
        method_outputs(METHODS, SAMPLE_IDS, config["output"])
    output:
        os.path.join(config["output"], "consensus.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(LOG_DIR,"consensus.log")
    params:
        methods = " ".join(METHODS)
    resources:
        mem_mb=50000,
        cpus_per_task=1
    shell:
        "python scripts/combine.py --input {config[output]} --methods {params.methods} --output {output} > {log} 2>&1"

rule report:
    input:
        consensus = os.path.join(config["output"], "consensus.csv"),
        samples = INPUT_FILES
    output:
        os.path.join(config["output"], "report.html")
    conda:
        "environment.yml"
    log:
        os.path.join(LOG_DIR,"report.log")
    resources:
        mem_mb=50000,
        cpus_per_task=1
    params:
        embedding = config.get("embedding", "spatial")
    shell:
        "python scripts/report.py --consensus {input.consensus} --input {config[input]} --output {output} --embedding {params.embedding} > {log} 2>&1"
