# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]


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
        """squigualiser reform --sig_move_offset {params.offset} --rna --kmer_length 1 -c --bam {input} -o {output}"""

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
        "squigualiser realign --rna --bam {input.bam} --paf {input.paf} -o {output}"

rule pod2slow:
    input:
        pod5 = "resources/pod5/{bam_file}/{pod5_file}.pod5"
    output:
        "resources/blow5/{bam_file}/{pod5_file}.blow5"
    conda:
        "../envs/bluecrab.yaml"
    threads: 16
    shell:
        "blue-crab p2s {input.pod5} -o {output}"

# Capture all blow5 files for the samples listed in samples.txt
files = expand("resources/blow5/{bam_file}/*.blow5", bam_file=samples)

rule signal2ref:
    input:
        slow5 = "resources/blow5/{bam_file}/{pod5_file}.blow5",
        realigned = "resources/alignments/squigualiser/{bam_file}_realigned_sorted.bam",
        realigned_index = "resources/alignments/squigualiser/{bam_file}_realigned_sorted.bam.bai",
        ref = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        temp = temp(directory("resources/.temp/{read_id}-{pod5_file}/{bam_file}_{region}")),
        html = "resources/signal/squigualiser/{read_id}-{pod5_file}/{bam_file}_{region}.html"
    conda:
        "../envs/squigualiser.yaml"
    params:
        chr = r"gi\|1154491913\|ref\|NR_003286.4\|"
    threads: 16
    shell:
        """
        squigualiser plot -r {wildcards.read_id} --rna --region {params.chr}:{wildcards.region} \
        --file {input.ref} --slow5 {input.slow5} --alignment {input.realigned} --output_dir {output.temp}; \
        mv {output.temp}/{wildcards.read_id}_.html {output.html}
        """
