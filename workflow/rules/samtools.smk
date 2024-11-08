# Load sample names from samples.txt
with open("samples.txt") as f:
    samples = [line.strip() for line in f if line.strip()]

rule samtools_sort:
    input:
        "resources/{dir}/{sample}.bam"
    output:
        "resources/{dir}/{sample}_sorted.bam"
    conda:  
        "../envs/samtools.yaml"
    threads: 16
    shell:
        "samtools sort -T resources/{wildcards.dir}/{wildcards.sample}_sorted -@ {threads} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "resources/{dir}/{sample}_sorted.bam"
    output:
        "resources/{dir}/{sample}_sorted.bam.bai"
    conda:  
        "../envs/samtools.yaml"
    threads: 16
    shell:
        "samtools index -@ {threads} {input}"

rule samtools_merge:
    input:
        expand("resources/{dir}/{sample}_sorted.bam", sample=samples, dir="basecalls")
    output:
        "resources/basecalls/wt_basecalls.bam"
    conda:
        "../envs/samtools.yaml"
    shell: 
        "samtools merge -o {output} {input}"
