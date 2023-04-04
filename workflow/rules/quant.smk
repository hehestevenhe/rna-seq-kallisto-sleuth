def genomebam_inputs(wildcards):
    activated = config["genomebam"]["activate"]
    if activated == True:
        rule_input = "resources/chrom_edit.txt"
    elif activated == False:
        rule_input = "dummy.txt"
    return rule_input
        
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
        
rule star_align:
    input:
        unpack(get_trimmed_star),  
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
    threads: 12
    resources: 
        tmpdir="./"
    wrapper:
        "https://github.com/hehestevenhe/rna-seq-kallisto-sleuth/raw/main/workflow/wrapper"
        
rule star_bam_naming:
    input:
        file="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        bam="results/star/indexed/{sample}-{unit}.bam",
    shell:
        "mv {input.file} results/star/indexed/{wildcards.sample}-{wildcards.unit}.bam"
        
rule star_bam_indexing:
    input: 
        "results/star/indexed/{sample}-{unit}.bam",
    output:
        bai="results/star/indexed/{sample}-{unit}.bam.bai",
    log: 
        "logs/star/{sample}-{unit}-indexing.log",
    wrapper:
        "v1.25.0/bio/samtools/index"

