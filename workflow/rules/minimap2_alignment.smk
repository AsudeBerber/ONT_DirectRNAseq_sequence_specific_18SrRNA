rule minimap2_align_txome:
    input:
        bam = "resources/basecalls/{sample}_basecalls.bam",
        fa = "resources/referencetranscriptome/18SrRNA.fa"
    output:
        "resources/alignments/{sample}_aligned.bam"
    conda:
        "../envs/minimap2.yaml"
    threads:
        16
    shell:
    # MM and ML include modified bases/probability, s. https://samtools.github.io/hts-specs/SAMtags.pdf, p.7
    # so far, this is not used further, but could be included with a function similar to access_mv() in align_signal.py
        """
        samtools fastq -@ {threads} -T mv,ts,ns,MM,ML {input.bam} |
            minimap2 -ax map-ont -k14 --secondary=no -t {threads} {input.fa} - -y --MD|
            samtools view -F 2048 -bh -@ {threads} -o {output}
        """ 

rule minimap2_align_txome_all:
    input:
        "resources/alignments/p2s_aligned.bam",
        "resources/alignments/p2i_aligned.bam",
        "resources/alignments/ko_aligned.bam"


## rule ....all input file and output file extensions should match with each other