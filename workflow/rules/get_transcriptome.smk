rule get_transcriptome:
    output:
        "resources/referencetranscriptome/18SrRNA.fa"
    params:
        url = "https://www.ncbi.nlm.nih.gov/nuccore/NR_003286.4?report=fasta" 
    shell:
        """
        wget -O {output} {params.url}
        """


# to get transcriptome file human whole transcriptome fasta file = use GENCODE
# because gencode transcriptome fasta file contains "chr1" "chr2".. ENSEMBL not and this can create a problem. 
# to make it reproducible always use the link of datafiles getting from internet. Do not download and use data or seq file from website.
