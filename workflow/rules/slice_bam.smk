# creates bam file with single read for testing purposes: 1. gets data for single read 2. samtools needs header -> concatenates header + read data
rule bam_single_read:
    input: "resources/alignments/p2s_aligned.bam"
    output: "resources/alignments/{read_ID}.bam"
    conda: "../envs/slice_bam.yaml"
    threads: 8
    shell:
        """samtools view {input} | grep {wildcards.read_ID} > resources/.temp/{wildcards.read_ID} \n
        samtools view -h {input} | head -n2 | cat - resources/.temp/{wildcards.read_ID} | samtools view -o {output} \n
        rm resources/.temp/{wildcards.read_ID}"""
