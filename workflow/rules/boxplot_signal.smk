EVENTS = list(range(0,9))
#dir without .npz ending
#dir = "CCG_window_21_p2s_aligned"
rule plot_boxplot:
    input: 
        "resources/results/p2s/{dir}.npz"
    output:
        svg = ["resources/signal/p2s/signal_summary/{dir}/1337_1842_430_event_{event}.svg".format(event = event) for event in EVENTS]
    conda:
        "../envs/boxplot.yaml"
    params:
        # how many bases to plot left/right to the middle base
        window = 8
    threads: 8
    shell:
        "python workflow/scripts/boxplot_signal.py -f {input} -w {params.window}"