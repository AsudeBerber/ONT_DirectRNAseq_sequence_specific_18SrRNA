# creates bam file with single read for testing purposes
rule bam_single_read:
    input: "resources/alignments/p2s_aligned.bam"
    output: "resources/alignments/{read_ID}.bam"
    conda: "../envs/slice_bam.yaml"
    threads: 8
    shell:
        "samtools view -h {input} | grep {wildcards.read_ID} | samtools view -h -o {output}"
