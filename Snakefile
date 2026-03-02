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
os.makedirs(os.path.join(config["output"], "logs"), exist_ok=True)

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
        os.path.join(config["output"], "logs", "singler_{sample}.log")
    resources:
        mem_mb=120000,
        cpus_per_task=20
    shell:
        "Rscript scripts/run_singler.R --input {input} --ref {config[ref]} --output {output} > {log} 2>&1"

rule tacco:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "tacco", "{sample}_tacco.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "tacco_{sample}.log")
    resources:
        mem_mb=50000,
        cpus_per_task=4
    shell:
        "python scripts/run_tacco.py --input {input} --ref {config[ref]} --output {output} > {log} 2>&1"

rule phispace:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "phispace", "{sample}_phispace.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "phispace_{sample}.log")
    resources:
        mem_mb=150000,
        cpus_per_task=4
    shell:
        "Rscript scripts/run_phi_space.R --input {input} --ref {config[ref]} --output {output} > {log} 2>&1"

rule rctd:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "rctd", "{sample}_rctd.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "rctd_{sample}.log")
    threads: 5
    resources:
        mem_mb=150000,
        cpus_per_task=5
    shell:
        "Rscript scripts/run_rctd.R --input {input} --ref {config[ref]} --output {output} --max_cores {threads} > {log} 2>&1"

rule insitutype:
    input:
        get_sample_input
    output:
        os.path.join(config["output"], "insitutype", "{sample}_insitutype.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "insitutype_{sample}.log")
    resources:
        mem_mb=50000,
        cpus_per_task=4
    shell:
        "Rscript scripts/run_insitutype.R --input {input} --ref {config[ref]} --output {output} > {log} 2>&1"

rule consensus:
    input:
        method_outputs(METHODS, SAMPLE_IDS, config["output"])
    output:
        os.path.join(config["output"], "consensus.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "consensus.log")
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
        os.path.join(config["output"], "logs", "report.log")
    resources:
        mem_mb=50000,
        cpus_per_task=1
    shell:
        "python scripts/report.py --consensus {input.consensus} --input {config[input]} --output {output} --embedding {config.get('embedding', 'spatial')} > {log} 2>&1"
