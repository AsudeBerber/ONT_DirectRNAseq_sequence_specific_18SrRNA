rule samtools_sort:
    input:
        "resources/{dir}/{sample}.bam"
    output:
        "resources/{dir}/{sample}_sorted.bam"
    conda:  
        "../envs/samtools.yaml"
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
        "../envs/samtools.yaml"
    threads:
        16
    shell:
        "samtools index -@ {threads} {input}"


rule samtools_merge:
    input:
        "resources/basecalls/p2i_basecalls_sorted.bam",
        "resources/basecalls/p2s_basecalls_sorted.bam"
    output:
        "resources/basecalls/wt_basecalls.bam"
    conda:
        "../envs/samtools.yaml"
    shell: 
        "samtools merge -o {output} {input}"
