rule reformat:
    input:
        "resources/alignments/{bam_file}.bam"
    output:
        "resources/alignments/squigualiser/{bam_file}_reform.paf"
    conda:
        "../envs/squigualiser.yaml"
    params:
        offset = 0
    threads: 16
    shell:
        """ squigualiser reform --sig_move_offset {params.offset} --rna --kmer_length 1 -c --bam {input} -o {output}"""

rule realign:
    input:
        paf = "resources/alignments/squigualiser/{bam_file}_reform.paf",
        bam = "resources/alignments/{bam_file}.bam"
    output:
        "resources/alignments/squigualiser/{bam_file}_realigned.bam"
    conda:
        "../envs/squigualiser.yaml"
    threads: 16
    shell:
        "squigualiser realign --rna --bam {input.bam}  --paf {input.paf} -o {output}"

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

files = glob_wildcards("resources/blow5/p2s/{file}.blow5")

rule signal2ref:
    input:
        slow5 = "resources/blow5/p2s/PAW35875_9fd38647_68d05f77_211.blow5", 
        realigned = "resources/alignments/squigualiser/{bam_file}_realigned.bam",
        realigned_index = "resources/alignments/squigualiser/{bam_file}_realigned.bam.bai",
        ref = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        "{output_dir}/{read_id}/{bam_file}_{region}.html"
    conda:
        "../envs/squigualiser.yaml"
    params:
        read_id = "7be77036-bcfb-4493-a95a-dd58b6975e5b",
        OUTPUT_DIR = "resources/signal/p2s/squigualiser/{params.read_id}",
        chr = r"gi\|1154491913\|ref\|NR_003286.4\|",
    threads: 16
    shell:
        """squigualiser plot  --rna --region {params.chr}:{wildcards.region}\
        --file {input.ref} --slow5 {input.slow5} --alignment {input.realigned} --output_dir resources/.temp \n 
        mv resources/.temp/{read_id}_.html {output}"""
