# creates bam file with single read for testing purposes: 1. gets data for single read 2. samtools needs header -> concatenates header + read data
rule bam_single_read:
    input: "resources/alignments/p2s_aligned.bam"
    output: 
        bam_temp = ("resources/.temp/{read_ID}.bam"),
        bam = "resources/alignments/{read_ID}.bam"
    conda: "../envs/slice_bam.yaml"
    threads: 8
    shell:
        """mkdir -p resources/.temp \n
        samtools view {input} | grep {wildcards.read_ID} > {output.bam_temp} \n
        samtools view -h {input} | head -n2 | cat - {output.bam_temp} | samtools view -o {output.bam}"""
