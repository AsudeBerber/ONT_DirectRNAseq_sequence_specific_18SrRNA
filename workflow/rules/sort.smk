rule samtools_sort:
    input:
        "resources/{dir}/{sample}.bam"
    conda:  
        "../envs/sort.yaml"
    output:
        "resources/{dir}/{sample}_sorted.bam"
    threads:
        16
    shell:
        "samtools sort -T data/mapped/{wildcards.dir}/{wildcards.sample}_sorted -@ {threads} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "resources/{dir}/{sample}.bam"
    conda:  
        "../envs/sort.yaml"
    output:
        "resources/{dir}/{sample}.bam.bai"
    threads:
        16
    shell:
        "samtools index -@ {threads} {input}"



