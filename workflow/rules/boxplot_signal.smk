rule plot_boxplot:
    input: 
        "resources/results/p2s/CCG_window_21_p2s_aligned_subsample_0001.npz"
    output:
        f"resources/signal/p2s/signal_summary/{dir_save}/1337_1842_430_event_{event}.svg"
    conda:
        "../envs/"
    params:
        # how many bases to plot left/right to the middle base
        window = 8
    shell:
        "python workflow/scripts/boxplot_signal.py -f {input} -w {params.window}"