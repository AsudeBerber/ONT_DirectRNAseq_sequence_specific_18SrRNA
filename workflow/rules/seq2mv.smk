# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

# Load read IDs from read_ids.txt
with open("read_ids.txt") as f:
    read_ids = [line.strip() for line in f if line.strip()]

# Define the position and range of interest (this can be adapted or read from a config if needed)
pos = 1337  # Example position
range_val = 100  # Example range value

rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/{sequencer}_aligned_sorted.bam",
        bai = "resources/alignments/{sequencer}_aligned_sorted.bam.bai"
    output:
        "resources/signal/{sequencer}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg"
    params: 
        region = r"gi\|1154491913\|ref\|NR_003286.4\|"
    wildcard_constraints:
        sequencer = "|".join(sequencers)
    conda:
        "../envs/seq2mv.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/seq2mv_direct_RNA.py \
            --sequencer {wildcards.sequencer} \
            --sample {input.bam} \
            --readID {wildcards.read_id} \
            --pos {wildcards.pos} --range {wildcards.range} \
            --pod5-dir resources/pod5/{wildcards.sequencer} \
            --region {params.region}
        """

rule seq2mv_single_read_all:
    input:
        expand(
            "resources/signal/{sequencer}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg",
            sequencer=sequencers,
            read_id=read_ids,
            pos=pos,
            range=range_val
        )
