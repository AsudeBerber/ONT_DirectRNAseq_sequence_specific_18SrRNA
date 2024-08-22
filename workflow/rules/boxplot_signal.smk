"""
boxplot_signal.py makes boxplots for 9 features given in .npz file (generated
by signal_summary.py), currently for {window} bases around positions 430 (no 
acetylation), 1337 and 1842 takes about ~1h to run on PromethIon, max. RAM usage
< 50%; multiple cores are used by Pandas/Matplotlib(?) motif: "CCG", not
relevant for code; window: total size of positions covered e.g. 21 (middle pos.
+ 10 left and right), event: events 0-8 (s. signal_summary.smk)
"""

EVENTS = list(range(0,9))

rule plot_boxplot:
    input: 
        "resources/signal_summary/CCG_window_{window_size,[0-9]+}_{sample}.npz"
    output:
        svg = expand("results/signal_summary/CCG_window_{{window_size,[0-9]+}}_{{sample}}/1337_1842_430_event_{event}.svg", event = EVENTS) 
    conda:
        "../envs/boxplot.yaml"
    threads: 
        16
    shell:
    # --no-mmap command disables mmap (loading to disk), could make problems on local PCs without enough RAM 
    # for some reason the mmap arg (https://numpy.org/doc/stable/reference/generated/numpy.memmap.html) has no effect when run on the promethion,
    # (takes ~40% RAM), but prevents crashing when computer memory is smaller than size of loaded file
        """
        set -x
        python workflow/scripts/boxplot_signal.py -f {input} -w {wildcards.window_size} --output-dir results/signal_summary/CCG_window_{wildcards.window_size}_{wildcards.sample}/
        
        """

rule plot_boxplot_all:
    input:
        expand("results/signal_summary/CCG_window_21_p2s/1337_1842_430_event_{event}.svg",
            event=EVENTS),
            
        expand("results/signal_summary/CCG_window_21_ko/1337_1842_430_event_{event}.svg",
            event=EVENTS),
        
        expand("results/signal_summary/CCG_window_21_p2i/1337_1842_430_event_{event}.svg",
            event=EVENTS)