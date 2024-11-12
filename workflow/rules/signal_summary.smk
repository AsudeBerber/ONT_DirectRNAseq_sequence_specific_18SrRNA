# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]


rule make_pod5_index:
    input: 
        "resources/pod5/{sample}/"
    output: 
        "resources/signal_summary/{sample}_pod5_index.json"
    conda:
        "../envs/signal_sum.yaml"
    threads: 1
    shell:
        "python workflow/scripts/pod5_index.py --pod5 {input} -o {output}"

rule signal_sum:
    input: 
        pod5 = "resources/pod5/{sample}/",
        bam = "resources/alignments/{sample}_aligned.bam",
        json = "resources/signal_summary/{sample}_pod5_index.json"
    output: 
        "resources/signal_summary/{motif}_window_{window_size,[0-9]+}_{sample}.npz"
    conda:
        "../envs/signal_sum.yaml"
    threads:
        1
    shell:
        "python workflow/scripts/signal_summary.py \
            --json {input.json} \
            --bam {input.bam} \
            --window {wildcards.window_size} \
            --output {output}"

rule signal_sum_all:
    input:
        expand("resources/signal_summary/CCG_window_21_{sample}.npz", sample=samples)
