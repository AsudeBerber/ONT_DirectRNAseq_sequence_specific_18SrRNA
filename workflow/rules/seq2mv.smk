# Read sample names from a file specified by the config parameter
samples_file = config["samples_file"]

with open(samples_file) as f:
    samples = [line.strip() for line in f if line.strip()]

# Define the position and range of interest (this can be adapted or read from a config if needed)
pos = 1337  # Example position
range_val = 100  # Example range value

def get_all_read_ids():
    """Retrieve read IDs for all samples."""
    all_read_ids = {}
    for sample in samples:
        read_ids_file = f"resources/{sample}_read_ids.txt"
        with open(read_ids_file) as f:
            all_read_ids[sample] = [line.strip() for line in f if line.strip()]
    return all_read_ids

all_read_ids = get_all_read_ids()

def generate_output_files():
    """Generate all output file paths for the seq2mv_single_read_all rule."""
    output_files = []
    for sample in samples:
        for read_id in all_read_ids[sample]:
            output_files.append(
                f"resources/signal/{sample}/plots/{read_id}/{read_id}_{pos}-pm{range_val}.svg"
            )
    return output_files

output_files = generate_output_files()

rule seq2mv_single_read:
    input: 
        bam = "resources/alignments/{sample}_sorted.bam",
        bai = "resources/alignments/{sample}_sorted.bam.bai",
        #read_ids_file = "resources/{sample}_read_ids.txt"
    output:
        "resources/signal/{sample}/plots/{read_id}/{read_id}_{pos}-pm{range}.svg"
    params: 
        region = r"gi\|1154491913\|ref\|NR_003286.4\|"
    wildcard_constraints:
        sample = "|".join(samples),
        read_id = "|".join(set(sum(all_read_ids.values(), [])))
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
        output_files
