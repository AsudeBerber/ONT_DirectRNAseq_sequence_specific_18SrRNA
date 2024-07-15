# boxplot_signal.py makes boxplots for 9 features given in .npz file (generated by signal_summary.py), currently for {window} bases around positions 430 (no acetylation), 1337 and 1842
# takes about ~1h to run on PromethIon, max. RAM usage < 50%; multiple cores are used by Pandas/Matplotlib(?)
#motif: "CCG", not relevant for code; window: total size of positions covered e.g. 21 (middle pos. + 10 left and right), event: events 0-8 (s. signal_summary.smk)
EVENTS = list(range(0,9))

rule plot_boxplot:
    input: 
        "resources/results/p2s/{motif}_window_{window}_{bam_file}.npz"
    output:
        svg = expand("resources/signal/p2s/signal_summary/{{motif,[A-Za-z]+}}_window_{{window,[0-9]+}}_{{bam_file}}/1337_1842_430_event_{event}.svg", event = EVENTS) 
    conda:
        "../envs/boxplot.yaml"
    threads: 16
    shell:
    # --no-mmap command disables mmap (loading to disk), could make problems on local PCs without enough RAM 
    # for some reason the mmap arg (https://numpy.org/doc/stable/reference/generated/numpy.memmap.html) has no effect when run on the promethion,
    # (takes ~40% RAM), but prevents crashing when computer memory is smaller than size of loaded file
        "python workflow/scripts/boxplot_signal.py -f {input} -w {wildcards.window} --output-dir {wildcards.motif}_window_{wildcards.window}_{wildcards.bam_file}"
