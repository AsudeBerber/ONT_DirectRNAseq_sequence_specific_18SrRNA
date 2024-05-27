rule get_transcriptome:
    input:
        FTP.remote("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz")
    output:
        ""
    shell:
        "gunzip {inout} > {output}"

#to get transcriptome file human whole transcriptome fasta file = use GENCODE
#because gencode transcriptome fasta file contains "chr1" "chr2".. ENSEMBL not and this can create a problem. 
#to make it reproducible always use the link of datafiles getting from internet. Do not download and use data or seq file from website.