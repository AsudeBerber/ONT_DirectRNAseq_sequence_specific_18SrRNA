# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

# Define the position and range of interest (this can be adapted or read from a config if needed)
pos = 1337  # Example position
range_val = 100  # Example range value

rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/{sample}_aligned_sorted.bam",
        bai = "resources/alignments/{sample}_aligned_sorted.bam.bai",
        read_ids_file = "resources/{sample}_read_ids.txt"
    output:
        "resources/signal/{sample}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg"
    params: 
        region = r"gi\|1154491913\|ref\|NR_003286.4\|"
    wildcard_constraints:
        sample = "|".join(samples)
    conda:
        "../envs/seq2mv.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/seq2mv_direct_RNA.py \
            --sample {wildcards.sample} \
            --sample {input.bam} \
            --readID {wildcards.read_id} \
            --pos {wildcards.pos} --range {wildcards.range} \
            --pod5-dir resources/pod5/{wildcards.sample} \
            --region {params.region}
        """

rule seq2mv_single_read_all:
    input:
        expand(
            "resources/signal/{sample}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg",
            sample=samples,
            read_id=lambda wildcards: get_read_ids(wildcards.sample),
            pos=pos,
            range=range_val
        )

def get_read_ids(sample):
    """Retrieve read IDs from the appropriate file for a given sample."""
    read_ids_file = f"resources/{sample}_read_ids.txt"
    with open(read_ids_file) as f:
        return [line.strip() for line in f if line.strip()]
