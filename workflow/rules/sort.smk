rule samtools_sort:
    input:
        "data/mapped/{dir}/{sample}.bam"
    conda:  
        "../envs/sort.yaml"
    output:
        "data/mapped/{dir}/{sample}_sorted.bam"
    threads:
        16
    shell:
        "samtools sort -T data/mapped/{wildcards.dir}/{wildcards.sample}_sorted -@ {threads} "
        "-O bam {input} > {output}"


#transcriptome: should have _sorted_transcripts.bam ending
rule samtools_index:
    input:
        "data/mapped/{dir}/{sample}.bam"
    conda:  
        "../envs/sort.yaml"
    output:
        "data/mapped/{dir}/{sample}.bam.bai"
    threads:
        16
    shell:
        "samtools index -@ {threads} {input}"

rule samtools_sort_transcriptome:
    input:
        "data/mapped/transcriptome/{sample}_unsorted.bam"
    conda:  
        "../envs/sort.yaml"
    output:
        "data/mapped/transcriptome/{sample}_sorted_transcripts.bam"
    threads:
        16
    shell:
        "samtools sort -o {output} -@ {threads} data/mapped/transcriptome/{wildcards.sample}_unsorted.bam"


