# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

print(f"Samples: {samples}")
rule samtools_sort:
    '''print(f"Sort: {sample}")'''
    input:
        "resources/alignments/{sample}_aligned.bam"
    output:
        "resources/alignments/{sample}_sorted.bam"
    conda:  
        "../envs/samtools.yaml"
    threads: 16
    shell:
        "samtools sort -T resources/alignments/{wildcards.sample}_sorted -@ {threads} -O bam {input} > {output}"
        

rule samtools_index:
    input:
        "resources/alignments/{sample}_sorted.bam"
    output:
        "resources/alignments/{sample}_sorted.bam.bai"
    conda:  
        "../envs/samtools.yaml"
    threads: 16
    shell:
        "samtools index -@ {threads} {input}"

rule samtools_all:
    input:
        expand("resources/alignments/{sample}_sorted.bam.bai", sample=samples)


'''
rule samtools_merge:
    input:
        expand("resources/{dir}/{sample}_sorted.bam", sample=samples, dir="basecalls")
    output:
        "resources/basecalls/wt_basecalls.bam"
    conda:
        "../envs/samtools.yaml"
    shell: 
        "samtools merge -o {output} {input}"
'''
