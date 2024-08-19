"""
creates dict of all reads with corresponding pod5 file, this makes signal_summary.py much faster (no looping needed)
"""
rule make_pod5_index:
    input: 
        "resources/pod5/{sample}/"
    output: 
        "resources/pod5/index/{sample}/pod5_index.json"
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

FEATURES:
Normalises and collapses the signal based on the moves table. Outputs an array with the
following values for each called based:
4: log10 signal length
5: mean signal intensity
6: standard deviation of signal intensity
7: median signal intensity
8: median absolute deviation of signal intensity
0-3: mean signal intensity for each quartile
"""
# https://github.com/hiruna72/squigualiser/tree/main/docs contains many useful explanations on the movetable, pore models etc.
rule signal_sum:
    input: pod5 = "resources/pod5/{sample}/",
           bam = "resources/alignments/{sample}_basecalls.bam",
           json = "resources/pod5/index/{sample}/pod5_index.json"
    output: 
        "resources/signal_summary/{motif}_window_{window_size,[0-9]+}_{sample}.npz"
    conda:
        "../envs/signal_sum.yaml"
    threads:
        1
    shell:
        "python workflow/scripts/signal_summary.py --pod5 {input.pod5} --bam {input.bam} --window {wildcards.window_size} --output {output}"

