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
        pod5 = "resources/pod5/p2s/{pod5_file}.pod5" # test.bam: resources/pod5/p2s/PAW35875_9fd38647_68d05f77_211.pod5
    output:
        "resources/blow5/p2s/{pod5_file}.blow5"
    conda:
        "../envs/bluecrab.yaml"
    threads: 16
    shell:
        "blue-crab p2s  {input.pod5} -o {output}"

rule signal2ref:
    input:
        slow5 = "resources/blow5/p2s/*", 
        realigned = "resources/alignments/squigle/{bam_file}_realigned.bam",
        ref = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        "{params.OUTPUT_DIR}/{bamfile}_{params.region}.svg"
    conda:
        "../envs/squigle.yaml"
    params:
        OUTPUT_DIR = "resources/signal/p2s/squigle",
        region = r"gi\|1154491913\|ref\|NR_003286.4\|1330:1350",
        tag = 'optionA'
    threads: 16
    shell:
        """squigualiser plot --file ${input.ref} --slow5 ${input.slow5} --alignment ${input.realigned} --output_dir ${params.OUTPUT_DIR} \
         --region ${params.region} --tag_name {params.tag} """
