# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]


rule bam_single_read:
    input: 
        bam = "resources/alignments/{sample}_aligned.bam"
    output: 
        bam_temp = temp("resources/.temp/{read_ID}.bam")
    conda: 
        "../envs/slice_bam.yaml"
    threads: 8
    shell:
        """set -x; mkdir -p resources/.temp resources/alignments/single_reads; \
        samtools view {input.bam} | grep {wildcards.read_ID} > {output.bam_temp}"""

rule bam_single_read_2:
    input: 
        bam = "resources/alignments/{sample}_aligned.bam",
        bam_temp = "resources/.temp/{read_ID}.bam"
    output: 
        bam = "resources/alignments/single_reads/{read_ID}.bam"
    conda: 
        "../envs/slice_bam.yaml"
    threads: 8
    shell:
        # Include header in the single-read BAM file
        """samtools view -h {input.bam} | head -n4 | cat - {input.bam_temp} | samtools view -bh -o {output} || true"""
