rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/{sequencer}_aligned_sorted.bam",
        bai = "resources/alignments/{sequencer}_aligned_sorted.bam.bai"
    output:
        "resources/signal/{sequencer}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg"
    params:
        read = "1fee0116-fdcc-4647-af43-9ea8d074de19", 
        # site 1: 1337	ac4C 79%
        # site 2: 1842	ac4C 99%
        # region = reference span (in IGV read description)
        # thousand seperators have to be removed (e.g. 1.868 -> 1868); special characters like "|" have to be written with escape sign ("|" -> "\|")
        region = r"gi|1154491913|ref|NR_003286.4|:15-1868"
    wildcard_constraints:
        sequencer = "p2i|p2s"
    conda:
        "../envs/seq2mv.yaml"
    shell:
     """python workflow/scripts/seq2mv_direct_RNA.py \
        {wildcards.sequencer} \
        {input.bam} \
        {wildcards.read_id} \
        {wildcards.pos} {wildcards.range} \
        --pod5_dir resources/pod5/{wildcards.sequencer} \
        --region {params.region}"""
     
 