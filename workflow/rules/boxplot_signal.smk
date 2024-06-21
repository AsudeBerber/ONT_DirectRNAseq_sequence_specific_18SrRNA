EVENTS = range(0,9)
rule plot_boxplot:
    input: 
        "resources/results/p2s/{dir}.npz"
    output:
        expand(f"resources/signal/p2s/signal_summary/{dir}/1337_1842_430_event_{event}.svg", event = EVENTS)
    conda:
        "../envs/"
    params:
        # how many bases to plot left/right to the middle base
        window = 8
    threads: 8
    shell:
        "python workflow/scripts/boxplot_signal.py -f {input} -w {params.window}"