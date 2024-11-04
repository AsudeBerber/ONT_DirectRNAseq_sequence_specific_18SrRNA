rule dorado:
    input:
        directory("resources/pod5/{sample}")
    output:
        "resources/basecalls/{sample}_basecalls.bam"
    params:
        model = "sup@v5.0.0",
        modification = ",pseU@v1,m6A@v1",
        batch_size = 320
    wildcard_constraints:
        sample = "ko_thn|p2i_thn"
    threads:
        20
    shell:
        "dorado_0.7.3 basecaller {params.model}{params.modification} -b {params.batch_size} {input} --emit-moves > {output}"


rule dorado_all:
    input:
        "resources/basecalls/p2i_basecalls.bam.thn.bam",
        "resources/basecalls/ko_basecalls.bam.thn.bam"



# merge these data p2i and p2s


    
# dorado_all apply for all rules since I have two samples data files
# while running it snakemake -np dorado_all (it recognizes all input files and rule dorado as well)
# then snakemake dorado_all --cores 4 (to run it actual)
