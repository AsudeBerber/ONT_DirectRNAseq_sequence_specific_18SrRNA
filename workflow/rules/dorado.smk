# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]


rule dorado:
    input:
        directory("resources/pod5/{sample}")
    output:
        "resources/basecalls/{sample}_basecalls.bam"
    params:
        model = "sup@v5.0.0",
        modification = ",pseU@v1,m6A@v1",
        batch_size = 320
    threads:
        20
    shell:
        "dorado_0.7.3 basecaller {params.model}{params.modification} -b {params.batch_size} {input} --emit-moves > {output}"

rule dorado_all:
    input:
        expand("resources/basecalls/{sample}_basecalls.bam", sample=samples)
