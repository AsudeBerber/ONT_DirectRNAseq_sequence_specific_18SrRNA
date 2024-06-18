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

directories, files = glob_wildcards("resources/blow5/p2s/{file}.pod5")

rule signal2ref:
    input:
        slow5 = expand("resources/blow5/p2s/{pod5_file}.pod5", pod5_file=files), 
        realigned = "resources/alignments/squigle/{bam_file}_realigned.bam",
        ref = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        "{output_dir}/{bam_file}_{start}_{stop}.svg"
    conda:
        "../envs/squigle.yaml"
    params:
        OUTPUT_DIR = "resources/signal/p2s/squigle",
        chr = r"gi\|1154491913\|ref\|NR_003286.4\|",
        region = "{wildcards.start}:{wildcards.stop}",
        tag = 'optionA'
    threads: 16
    shell:
        """squigualiser plot --file ${input.ref} --slow5 ${input.slow5} --alignment ${input.realigned} --output_dir ${wildcards.output_dir} \
         --region ${params.chr}{params.region} --tag_name {params.tag} """
