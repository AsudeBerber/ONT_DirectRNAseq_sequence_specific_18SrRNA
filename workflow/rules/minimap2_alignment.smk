rule minimap2_align_txome:
    input:
        bam = "resources/basecalls/{sample}_basecalls.bam",
        fa = "resources/referencetranscriptome/gencode.v46.transcripts.fa"
    output:
        "resources/alignments/{sample}_aligned_txome.bam"
    conda:
        "../envs/minimap2.yaml"
    threads:
        4
    shell:
        """
        samtools fastq -@ {threads} -T mv,ts,ns {input.bam} |
            minimap2 -ax map-ont -N 10 -k14 -t {threads} {input.fa} - |
            samtools view -bh -@ {threads} -o {output}
        """ 

rule minimap2_align_txome_all:
    input:
        "resources/basecalls/p2i_basecalls.bam",
        "resources/basecalls/p2s_basecalls.bam"