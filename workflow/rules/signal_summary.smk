"""
creates dict of all reads with corresponding pod5 file, this makes signal_summary.py much faster (no looping needed)
"""
rule make_pod5_index:
    input: 
        "resources/pod5/p2s/"
    output: 
        "resources/pod5/index/p2s/pod5_index.json"
    conda:
        "../envs/json.yaml"
    threads: 1
    shell:
        "python workflow/scripts/pod5_index.py --pod5 {input} -o pod5_index"


"""
calculates 9 features (s. align_signal.py) for given window_size (21 works very well, for 9-mers at least 17) for all reads in bam file at pos 429 (no acetylation, rel. similar sequence),
1337 and 1842 (this can be changed in line 27 of signal_summary.py)
! filtering of the bam_file for only the wanted reference transcript could be needed when aligning to several reference transcripts (currently the function only looks for the reference position,
not the "chromosome")
as this function calls align_signal.py, both functions are required to be in the same directory
"""
rule signal_sum:
    input: pod5 = "resources/pod5/p2s/",
           bam = "resources/alignments/{bam_file}.bam",
           json = "resources/pod5/index/p2s/pod5_index.json"
    output: 
        "resources/results/p2s/{motif}_window_{window_size,[0-9]+}_{bam_file}.npz"
    conda:
        "../envs/signal_sum.yaml"
    threads:
        1
    shell:
        "python workflow/scripts/signal_summary.py --pod5 {input.pod5} --bam {input.bam} --window {wildcards.window_size} --output {output}"

