rule dorado:
    input:
        directory("resources/pod5/{sample}")
    output:
        bam = "resources/basecalls/{sample}_basecalls.bam"
    params:
        model = "sup@v3.0.1"
    threads:
        4
    shell:
        "dorado_0.7.0 basecaller {params.model} {input} --emit-moves > {output.bam}"


rule dorado_all:
    input:
        ["resources/basecalls/p2i_basecalls.bam",
            "resources/basecalls/p2s_basecalls.bam"]
