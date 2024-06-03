rule minimap2_align_txome:
    input:
        bam = "resources/basecalls/{sample}_basecalls.bam",
        fa = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        "resources/alignments/{sample}_aligned.bam"
    conda:
        "../envs/minimap2.yaml"
    threads:
        8
    shell:
        """
        samtools fastq {input.bam} -T mv,ts,ns -@ {threads} \
        | minimap2 -ax map-ont -k14 --secondary=no -t {threads} {input.fa} - \
        | samtools view -F 2048 -bh -@ {threads} -o {output}
        """ 

rule minimap2_align_txome_all:
    input:
        "resources/alignments/p2i_aligned.bam",
        "resources/alignments/p2s_aligned.bam"
