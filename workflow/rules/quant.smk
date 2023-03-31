def genomebam_inputs(wildcards):
    activated = config["genomebam"]["activate"]
    if (activated == lower("true")):
        input = "resources/chrom_edit.txt"
    else:
        input = "dummy.txt"
    return input
        
ruleorder: kallisto_genomebam > kallisto_quant

rule kallisto_index:
    input:
        "resources/transcriptome.cdna.fasta",
    output:
        "results/kallisto/transcripts.idx",
    log:
        "results/logs/kallisto/index.log",
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} 2> {log}"


rule kallisto_quant:
    input:
        fq=get_trimmed,
        idx="results/kallisto/transcripts.idx",
    output:
        directory("results/kallisto/{sample}-{unit}"),
    log:
        "results/logs/kallisto/quant/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output} "
        "{params.extra} {input.fq} 2> {log}"

rule kallisto_genomebam:
    input:
        fq=get_trimmed,
        idx="results/kallisto/transcripts.idx",
        gtf="resources/genome.gtf",
        chrom=genomebam_inputs,
    output:
        dir=directory("results/kallisto/{sample}-{unit}"),
        bam="results/kallisto/{sample}-{unit}/pseudoalignments.bam",
    log:
        "results/logs/kallisto/genomebam/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    resources:
        tmpdir="./"
    shell:
        "kallisto quant -i {input.idx} -o {output.dir} "
        "{params.extra} --genomebam --gtf {input.gtf} --chromosomes {input.chrom} {input.fq} 2> {log}"
