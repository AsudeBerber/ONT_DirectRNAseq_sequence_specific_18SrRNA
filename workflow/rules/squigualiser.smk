rule reformat:
    input:
        "resources/alignments/{bam_file}.bam"
    output:
        "resources/alignments/squigle/{bam_file}_reform.paf"
    conda:
        "../envs/squigle.yaml"
    params:
        offset = 0
    threads: 16
    shell:
        """ squigualiser reform --sig_move_offset {params.offset} --kmer_length 1 -c --bam {input} -o {output}"""

rule realign:
    input:
        paf = "resources/alignments/squigle/{bam_file}_reform.paf",
        bam = "resources/alignments/{bam_file}.bam"
    output:
        "resources/alignments/squigle/{bam_file}_realigned.bam"
    conda:
        "../envs/squigle.yaml"
    threads: 16
    shell:
        "squigualiser realign --bam {input.bam} --paf {input.paf} -o {output}"

rule pod2slow:
    input:
        pod5 = "resources/pod5/p2s/{pod5_file}.pod5"
    output:
        "resources/blow5/p2s/{pod5_file}.blow5"
    conda:
        "../envs/bluecrab.yaml"
    threads: 16
    shell:
        "blue-crab p2s {pod5.input} -o {output}"

rule signal2ref:
    input:
        slow5 = "", 
        realigned = "resources/alignments/squigle/{bam_file}_realigned.bam"
    output:
        "{OUTPUT_DIR}/output"
    conda:
        "../envs/squigle.yaml"
    threads: 16
    shell:
        """squigualiser plot --file ${REF} --slow5 ${input.slow5} --alignment ${input.realigned} --output_dir ${OUTPUT_DIR} \
         --region ${REGION} --tag_name 'optionA' """

"