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

# as this is only intended to do for a few reads, i would advise to check the corresponding pod5 file to the read_id (to find in resources/results/p2s/pod5.json)
# it would be possible to create one big slow5 file, however this would take up ~100 GB of storage
rule signal2ref:
    input:
        slow5 = "resources/blow5/p2s/{pod5_file}.blow5", 
        realigned = "resources/alignments/squigualiser/{bam_file}_realigned_sorted.bam",
        realigned_index = "resources/alignments/squigualiser/{bam_file}_realigned_sorted.bam.bai",
        ref = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        # e.g resources/signal/squigualizer/READ_ID/p2s_aligned_1800-1850.html
        temp = directory("$TEMP/squigualiser/{read_id}-{pod5_file}/{bam_file}_{region}"),
        html = "resources/signal/squigualiser/{read_id}-{pod5_file}/{bam_file}_{region}.html"
    conda:
        "../envs/squigualiser.yaml"
    params:
        # this is the "chromosome" for 18S rRNA, change otherwise
        chr = r"gi\|1154491913\|ref\|NR_003286.4\|"
    threads: 16
    shell:
        """mkdir -p {output.temp}; squigualiser plot -r {wildcards.read_id} --rna --region {params.chr}:{wildcards.region}\
        --file {input.ref} --slow5 {input.slow5} --alignment {input.realigned} --output_dir {output.temp}; \
        mv {output.temp}/{wildcards.read_id}_.html {output.html}; rm -r {output.temp}"""
