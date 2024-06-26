rule signal_sum:
    input: pod5 = "resources/pod5/p2s/",
           bam = "resources/alignments/{bam_file}.bam"
    output: 
        "resources/results/p2s/{motif}_window_{window_size}_{bam_file}.npz"
    conda:
        ""../envs/signal_sum.yaml"
    #params:
        # window_size = 21
        # motif = "CCG"
    threads:
        16
    shell:
        "signal_summary.py --pod5 {params.pod5} --bam {input.bam} --window{params.window_size} --motif {params.motif} --output {output}"

