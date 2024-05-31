rule get_transcriptome:
    output:
        "resources/referencetranscriptome/18SrRNA.fa"
    params:
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1154491913&extrafeat=null&conwithfeat=on&hide-cdd=on&showgi=1&ncbi_phid=CE8B3062659506510000000008F20712"
    shell:
        """
        wget -O '{output}' {params.url}
        """


# to get transcriptome file human whole transcriptome fasta file = use GENCODE
# because gencode transcriptome fasta file contains "chr1" "chr2".. ENSEMBL not and this can create a problem. 
# to make it reproducible always use the link of datafiles getting from internet. Do not download and use data or seq file from website.


