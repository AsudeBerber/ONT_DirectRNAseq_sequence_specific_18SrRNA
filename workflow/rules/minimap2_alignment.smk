# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]


rule minimap2_align_txome:
    input:
        bam = "resources/basecalls/{sample}_basecalls.bam",
        fa = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        "resources/alignments/{sample}_aligned.bam"
    conda:
        "../envs/minimap2.yaml"
    threads: 16
    shell:
        """
        samtools fastq -@ {threads} -T mv,ts,ns,MM,ML {input.bam} |
            minimap2 -ax map-ont -k14 --secondary=no -t {threads} {input.fa} - -y --MD |
            samtools view -F 2048 -bh -@ {threads} -o {output}
        """

rule minimap2_align_txome_all:
    input:
        expand("resources/alignments/{sample}_aligned.bam", sample=samples)
