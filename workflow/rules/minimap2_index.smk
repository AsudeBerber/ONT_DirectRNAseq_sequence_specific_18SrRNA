rule get_transcriptome:
    input:
        FTP.remote("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz")
    output:
        "resources/referencetranscriptome/gencode.v46.transcripts.fa"
    shell:
        "gunzip {input} > {output}"

#to get transcriptome file human whole transcriptome fasta file = use GENCODE
#because gencode transcriptome fasta file contains "chr1" "chr2".. ENSEMBL not and this can create a problem. 
#to make it reproducible always use the link of datafiles getting from internet. Do not download and use data or seq file from website.


rule minimap2_index:
    input:
        "resources/referencetranscriptome/gencode.v46.transcripts.fa"
    output:
        "resources/referencetranscriptome/transcriptome_minimap2.mmi"
    conda:
        "../envs/minimap2.yaml"
    threads:
        3
    shell:
        "minimap2 -ax splice -uf -k14 -t {threads} -d {output} {input}"