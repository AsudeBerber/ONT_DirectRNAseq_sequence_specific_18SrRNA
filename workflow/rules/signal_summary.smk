rule make_pod5_index:
    input: 
        "resources/pod5/p2s/"
    output: 
        "resources/pod5/index/p2s/pod5_index.json"
    conda:
        "../envs/json.yaml"
    threads: 8
    shell:
        "python workflow/scripts/pod5_index.py --pod5 {input} -o pod5_index"

rule signal_sum:
    input: pod5 = "resources/pod5/p2s/",
           bam = "resources/alignments/{bam_file}.bam",
           json = "resources/results/p2s/pod5_index.json"
    output: 
        "resources/results/p2s/{motif}_window_{window_size,[0-9]+}_{bam_file}.npz"
    conda:
        "../envs/signal_sum.yaml"
    #params:
        # window_size = 21
        # motif = "CCG"
    threads:
        16
    shell:
        "python workflow/scripts/signal_summary.py --pod5 {input.pod5} --bam {input.bam} --window {wildcards.window_size} --output {output}"

