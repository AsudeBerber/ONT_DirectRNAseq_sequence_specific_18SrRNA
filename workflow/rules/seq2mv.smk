rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/p2s_aligned_sorted.bam",
        bai = "resources/alignments/p2s_aligned_sorted.bam.bai"
    output:
        "resources/signal/{sequencer}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg"
    params:
        read = "7be77036-bcfb-4493-a95a-dd58b6975e5b", 
        # site 1: 1337	ac4C 79%
        # site 2: 1842	ac4C 99%
        # region = reference span (in IGV read description)
        # thousand seperators have to be removed (e.g. 1.868 -> 1868); special characters like "|" have to be written with escape sign ("|" -> "\|")
        region = r"gi\|1154491913\|ref\|NR_003286.4\|"
        # :15-1868"
    wildcard_constraints:
        sequencer = "p2i|p2s"
    conda:
        "../envs/seq2mv.yaml"
    shell:
     """python workflow/scripts/seq2mv_direct_RNA.py \
        --sequencer {wildcards.sequencer} \
        --sample {input.bam} \
        --readID {wildcards.read_id} \
        --pos {wildcards.pos} --range {wildcards.range} \
        --pod5_dir resources/pod5/{wildcards.sequencer} \
        --region {params.region}"""
     
 