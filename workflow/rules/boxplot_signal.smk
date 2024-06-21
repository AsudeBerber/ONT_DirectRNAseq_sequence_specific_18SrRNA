# takes about ~1h to run on PromethIon, max. RAM usage < 50%; multiple cores are used by Pandas/Matplotlib(?)
EVENTS = list(range(0,9))
#dir without .npz ending
rule plot_boxplot:
    input: 
        "resources/results/p2s/{dir}.npz"
    output:
        svg = expand("resources/signal/p2s/signal_summary/{dir}/1337_1842_430_event_{event}.svg", event = EVENTS, dir = "{dir}") 
    conda:
        "../envs/boxplot.yaml"
    params:
        # how many bases to plot left/right to the middle base
        window = 8
    threads: 16
    shell:
    # --no-mmap command disables mmap (loading to disk), could make problems on local PCs without enough RAM
        "python workflow/scripts/boxplot_signal.py -f {input} -w {params.window}"