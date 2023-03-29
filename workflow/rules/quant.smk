ruleorder: kallisto_genomebame > kallisto_quant

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
        chrom="resources/chrom_edit.txt",
    output:
        directory("results/kallisto/{sample}-{unit}"),
        directory("results/kallisto/{sample}-{unit}/pseudoalignment.bam"),
    log:
        "results/logs/kallisto/genomebam/{sample}-{unit}.log",
    params:
        extra=kallisto_params,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {output} "
        "{params.extra} --genomebam --gtf {input.gtf} --chromosomes {input.chrom} {input.fq} 2> {log}"
