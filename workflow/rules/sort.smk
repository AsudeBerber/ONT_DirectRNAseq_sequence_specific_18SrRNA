rule samtools_sort:
    input:
        "resources/{dir}/{sample}.bam"
    output:
        "resources/{dir}/{sample}_sorted.bam"
    conda:  
        "../envs/sort.yaml"
    threads:
        16
    shell:
        "samtools sort -T resources/{wildcards.dir}/{wildcards.sample}_sorted -@ {threads} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "resources/{dir}/{sample}.bam"
    output:
        "resources/{dir}/{sample}.bam.bai"
    conda:  
        "../envs/sort.yaml"
    threads:
        16
    shell:
        "samtools index -@ {threads} {input}"



