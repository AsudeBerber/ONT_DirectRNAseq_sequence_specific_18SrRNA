# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

EVENTS = list(range(0, 9))

rule plot_boxplot:
    input: 
        "resources/signal_summary/CCG_window_{window_size,[0-9]+}_{sample}.npz"
    output:
        svg = expand("results/signal_summary/CCG_window_{{window_size}}_{{sample}}/1337_1842_430_event_{event}.svg", event=EVENTS) 
    conda:
        "../envs/boxplot.yaml"
    threads: 
        16
    shell:
        # Command for generating boxplot SVGs
        "python workflow/scripts/boxplot_signal.py -f {input} -w {wildcards.window_size} --output-dir results/signal_summary/CCG_window_{wildcards.window_size}_{wildcards.sample}/"

rule plot_boxplot_all:
    input:
        expand("results/signal_summary/CCG_window_21_{sample}/1337_1842_430_event_{event}.svg", 
               sample=samples, event=EVENTS)
