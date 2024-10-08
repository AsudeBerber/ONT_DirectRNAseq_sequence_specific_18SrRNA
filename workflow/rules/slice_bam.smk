# creates bam file with single read for testing purposes: 1. gets data for single read 2. samtools needs header -> concatenates header + read data
rule bam_single_read:
    input: 
        "resources/alignments/p2s_aligned.bam"
    output: 
        bam_temp = temp("resources/.temp/{read_ID}.bam")
    conda: 
        "../envs/slice_bam.yaml"
    threads: 8
    shell:
        """set -x; mkdir -p resources/.temp resources/alignments/single_reads; \
        samtools view {input} | grep {wildcards.read_ID} > {output.bam_temp}"""

rule bam_single_read_2:
    input: 
        bam = "resources/alignments/p2s_aligned.bam",
        bam_temp = "resources/.temp/{read_ID}.bam"
    output: 
        bam = "resources/alignments/single_reads/{read_ID}.bam"
    conda: 
        "../envs/slice_bam.yaml"
    threads: 8
    shell:
        # bamfile needs header to be read by samtools, so the header is appended to the bam file of the single read
        # snakemake throws an error although this is working, ||true ignores all errors; this might be an issue with bash strict mode?
        """samtools view -h {input.bam} | head -n4 | cat - {input.bam_temp} | samtools view -bh -o {output} ||true"""
