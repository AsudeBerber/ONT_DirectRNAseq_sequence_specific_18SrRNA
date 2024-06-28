# boxplot_signal.py makes boxplots for 9 features given in .npz file (generated by signal_summary.py), currently for {window} bases around positions 430 (no acetylation), 1337 and 1842
# takes about ~1h to run on PromethIon, max. RAM usage < 50%; multiple cores are used by Pandas/Matplotlib(?)
EVENTS = list(range(0,9))

rule plot_boxplot:
    input: 
        "resources/results/p2s/{motif}_window_{window,[0-9]+}_{bam_file}.npz"
    output:
        svg = expand("resources/signal/p2s/signal_summary/{motif,[A-Za-z]}_window_{window,[0-9]+}_{bam_file}/1337_1842_430_event_{event}.svg", event = EVENTS, dir = "{dir}") 
    conda:
        "../envs/boxplot.yaml"
    threads: 16
    shell:
    # --no-mmap command disables mmap (loading to disk), could make problems on local PCs without enough RAM
        "python workflow/scripts/boxplot_signal.py -f {input} -w {params.window} --output-dir {wildcards.dir}"
