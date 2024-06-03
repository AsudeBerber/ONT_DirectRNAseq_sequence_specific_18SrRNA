rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/p2s_aligned_sorted.bam"
    output:
        "resources/signal/{sequencer}/plots/{read_id}/{read_id}_{start}-{end}.png"
    params:
        read = "1fee0116-fdcc-4647-af43-9ea8d074de19", 
        base_pos = "",
        bases_around =  "",
        # region = reference span (in IGV read description)
        region = "gi|1154491913|ref|NR_003286.4|:15-1.868"
    wildcard_constraints:
        sequencer = "p2i|p2s"
    conda:
        "../envs/seq2mv.yaml"
    shell:
     """python workflow/scripts/seq2mv_direct_RNA.py {input.bam} {wildcards.read_id} {wildcards.start} {wildcards.end} \
        --pod5_dir resources/pod5/{wildcards.sequencer} --region {params.region}"""
     
 