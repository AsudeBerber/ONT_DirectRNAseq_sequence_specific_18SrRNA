rule split_bam:
    input:
        "resources/{dir}/{bam_file}.bam"
    output:
        expand("{bam_file}_temp_{nr}.bam",
            nr=[1, 2, 3, 4, 5, 6, 7, 8])
    shell:
        "samtools view {input} | split -l 550000 > {output_dir}"

rule append_header:
    input:
      "split output"
    output:
        "{bam_file}_part_[1:8].bam"
    shell:
        "samtools view {input} -h | head -2 | cat - {wildcards.input} | samtools view -bh -o {wildcards.output}"