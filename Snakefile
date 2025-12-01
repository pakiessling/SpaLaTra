import glob
import os

INPUT_FILES = glob.glob(os.path.join(config["input"], '*.h5ad'))
SAMPLE_IDS = [os.path.splitext(os.path.basename(f))[0] for f in INPUT_FILES]


os.makedirs(os.path.join(config["output"],"phispace"), exist_ok=True)
os.makedirs(os.path.join(config["output"],"tacco"), exist_ok=True)
os.makedirs(os.path.join(config["output"],"singler"), exist_ok=True)
os.makedirs(os.path.join(config["output"],"rctd"), exist_ok=True)
os.makedirs(os.path.join(config["output"],"insitutype"), exist_ok=True)

localrules: all

rule all:
    input:
        expand(os.path.join(config["output"],"phispace", "{sample}_phispace.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"],"tacco", "{sample}_tacco.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "singler","{sample}_singler.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "rctd","{sample}_rctd.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "insitutype","{sample}_insitutype.csv"), sample=SAMPLE_IDS),
        os.path.join(config["output"],"consensus.csv")

    shell:
        """
        echo Annotation complete!
        """

rule singler:
    input:
        config["input"] + "/{sample}.h5ad"
    output:
        os.path.join(config["output"],"singler", "{sample}_singler.csv")
    conda:
        "enviroment.yml"
    shell:
        "Rscript scripts/run_singler.R --input {input} --ref {config[ref]} --output {output}"

rule tacco:
    input:
        config["input"] + "/{sample}.h5ad"
    output:
        os.path.join(config["output"],"tacco", "{sample}_tacco.csv")
    conda:
        "enviroment.yml"
    shell:
        "python scripts/run_tacco.py --input {input} --ref {config[ref]} --output {output}"

rule phispace:
    input:
        config["input"] + "/{sample}.h5ad"
    output:
        os.path.join(config["output"],"phispace", "{sample}_phispace.csv")
    conda:
        "enviroment.yml"
    shell:
        "Rscript scripts/run_phi_space.R --input {input} --ref {config[ref]} --output {output}"

rule rctd:
    input:
        config["input"] + "/{sample}.h5ad"
    output:
        os.path.join(config["output"],"rctd", "{sample}_rctd.csv")
    conda:
        "enviroment.yml"
    shell:
        "Rscript scripts/run_rctd.R --input {input} --ref {config[ref]} --output {output}"

rule insitutype:
    input:
        config["input"] + "/{sample}.h5ad"
    output:
        os.path.join(config["output"],"insitutype", "{sample}_insitutype.csv")
    conda:
        "enviroment.yml"
    shell:
        "Rscript scripts/run_insitutype.R --input {input} --ref {config[ref]} --output {output}"

rule consensus:
    input:
        expand(os.path.join(config["output"],"phispace", "{sample}_phispace.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"],"tacco", "{sample}_tacco.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "singler","{sample}_singler.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "rctd","{sample}_rctd.csv"), sample=SAMPLE_IDS),
        expand(os.path.join(config["output"], "insitutype","{sample}_insitutype.csv"), sample=SAMPLE_IDS)
    output:
        os.path.join(config["output"],"consensus.csv")
    conda:
        "enviroment.yml"
    shell:
        "python scripts/combine.py --input {config[output]} --output {output}"
