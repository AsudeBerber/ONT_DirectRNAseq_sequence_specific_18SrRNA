rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/{sample}_sorted.bam",
        bai = "resources/alignments/{sample}_sorted.bam.bai"
    output:
        plot = "resources/signal/ps2/plots/{read_id}/{read_id}_{pos}-pm{range}.svg",
        read_id_list = "resources/read_id_list_bam.txt"
    params: 
        region = r"gi\|1154491913\|ref\|NR_003286.4\|"
    conda:
        "../envs/seq2mv.yaml"
    threads: 1
    shell:
        """python workflow/scripts/seq2mv_direct_RNA.py \
           --sequencer {wildcards.sequencer} \
           --sample {input.bam} \
           --readID {wildcards.read_id} \
           --pos {wildcards.pos} --range {wildcards.range} \
           --pod5-dir resources/pod5/{wildcards.sequencer} \
           --region {params.region} \
           --output-read-list {output.read_id_list}"""

rule seq2mv_single_read_all:
    input:
        expand(
            "resources/signal/ps2/plots/{read_id}/{read_id}_{pos}-pm{range}.svg",
            read_id=get_read_ids(),
            pos=["1337", "1842"],
            range=["50"]
        )
