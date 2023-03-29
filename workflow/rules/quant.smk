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
        chrom="resources/chrom_edit.txt",
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
        
rule star_align:
    input:
        unpack(get_trimmed),  ### THIS WILL NEED SOME EXTRA WORK; NOT SURE IF MATCHES DESIRED INPUT FORMAT
        index="resources/star_genome",
    output:
        aln="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        idx=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "v1.21.4/bio/star/align"
        
