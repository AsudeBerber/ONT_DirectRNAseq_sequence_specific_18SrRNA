rule dorado:
    input:
        directory("resources/pod5/{sample}")
    output:
        "resources/basecalls/{sample}_basecalls.bam"
    params:
        model = "sup@v5.0.0",
        modification = ",m6A,pseU@v1"
    threads:
        20
    shell:
        "dorado_0.7.2 basecaller {params.model}{params.modification} {input} --emit-moves > {output}"


rule dorado_all:
    input:
        "resources/basecalls/p2i_basecalls.bam",
        "resources/basecalls/p2s_basecalls.bam"


#dorado_all apply for all rules since I have two samples data files
#while running it snakemake -np dorado_all (it recognizes all input files and rule dorado as well)
#then snakemake dorado_all --cores 4 (to run it actual)