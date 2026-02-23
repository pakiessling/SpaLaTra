import glob
import os

INPUT_FILES = glob.glob(os.path.join(config["input"], '*.h5ad'))
SAMPLE_IDS = [os.path.splitext(os.path.basename(f))[0] for f in INPUT_FILES]

os.makedirs(os.path.join(config["output"], "phispace"), exist_ok=True)
os.makedirs(os.path.join(config["output"], "tacco"), exist_ok=True)
os.makedirs(os.path.join(config["output"], "singler"), exist_ok=True)
os.makedirs(os.path.join(config["output"], "rctd"), exist_ok=True)
os.makedirs(os.path.join(config["output"], "insitutype"), exist_ok=True)
os.makedirs(os.path.join(config["output"], "logs"), exist_ok=True)

localrules: all

rule all:
    input:
        expand(os.path.join(config["output"], "phispace", "{sample}_phispace.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "tacco", "{sample}_tacco.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "singler", "{sample}_singler.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "rctd", "{sample}_rctd.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "insitutype", "{sample}_insitutype.csv"), sample=SAMPLE_IDS),
        os.path.join(config["output"], "consensus.csv"),
        os.path.join(config["output"], "report.html")
    shell:
        """
        echo Annotation complete!
        """

rule singler:
    input:
        config["input"] + "/{sample}.h5ad"
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
        config["input"] + "/{sample}.h5ad"
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
        config["input"] + "/{sample}.h5ad"
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
        config["input"] + "/{sample}.h5ad"
    output:
        os.path.join(config["output"], "rctd", "{sample}_rctd.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "rctd_{sample}.log")
    resources:
        mem_mb=150000,
        cpus_per_task=5
    threads: resources.cpus_per_task
    shell:
        "Rscript scripts/run_rctd.R --input {input} --ref {config[ref]} --output {output} --max_cores {threads} > {log} 2>&1"

rule insitutype:
    input:
        config["input"] + "/{sample}.h5ad"
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
        expand(os.path.join(config["output"], "phispace", "{sample}_phispace.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "tacco", "{sample}_tacco.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "singler", "{sample}_singler.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "rctd", "{sample}_rctd.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "insitutype", "{sample}_insitutype.csv"), sample=SAMPLE_IDS)
    output:
        os.path.join(config["output"], "consensus.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(config["output"], "logs", "consensus.log")
    resources:
        mem_mb=50000,
        cpus_per_task=1
    shell:
        "python scripts/combine.py --input {config[output]} --output {output} > {log} 2>&1"

rule report:
    input:
        consensus = os.path.join(config["output"], "consensus.csv"),
        samples = expand(config["input"] + "/{sample}.h5ad", sample=SAMPLE_IDS)
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
