rule reformat:
    input:
        "resources/alignments/{bam_file}.bam"
    output:
        "resources/alignments/PAF/{bam_file}_reform.paf"
    conda:
        "../envs/squigle.yaml"
    params:
        offset = 0
    threads: 16
    shell:
        """ squigualiser reform --sig_move_offset {params.offset} --kmer_length 1 -c --bam {input} -o {output}"""