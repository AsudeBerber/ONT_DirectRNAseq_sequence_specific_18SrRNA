rule dorado:
    input:
        directory("resources/pod5/{sample}")
    output:
        "resources/basecalls/{sample}_basecalls.bam"
    params:
        model = "sup@v3.0.1"
    threads:
        20
    shell:
        "dorado_0.7.0 basecaller {params.model} {input} --emit-moves > {output}"


rule dorado_all:
    input:
        "resources/basecalls/p2i_basecalls.bam",
        "resources/basecalls/p2s_basecalls.bam"


#dorado_all apply for all rules since I have two samples data files