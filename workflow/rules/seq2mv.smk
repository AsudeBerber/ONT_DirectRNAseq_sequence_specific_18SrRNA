"""
seq2mv_direct_RNA.py produces a plot of the sequencing signal (squigle) for a given read and around a given position in the reference +- range
sites of interest in 18S rRNA:     site 1: 1337	ac4C 79%        site 2: 1842	ac4C 99%
reads (and corresponding regions) can be found in the IGV viewer 
    or  (w/o region; gi... with no position stays the same for all 18S rRNA) in resources/pod5/index/p2s/pod5_index.json (can be generated via snakemake)


"""
rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/p2s_aligned_sorted.bam",
        bai = "resources/alignments/p2s_aligned_sorted.bam.bai"
    output:
        "resources/signal/{sequencer}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg"
    params: 
        # region = reference span (in IGV read description); e.g. "gi\|1154491913\|ref\|NR_003286.4\|:15-1868" (this would work only if the read exactly matches this, therefore it is also much faster)
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
        --pod5-dir resources/pod5/{wildcards.sequencer} \
        --region {params.region}"""
     
 