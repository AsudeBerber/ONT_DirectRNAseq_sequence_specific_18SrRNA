rule signal_sum_temp:
    input: pod5 = "",
           bam = ""
    output: 
        csv 1-7
    wildcards.args:
        pod5 parts
    conda:
        ""../envs/signal_sum.yaml"
    threads:
        16
    shell:
        "signal_summary.py"
rule merge_temp_files:
    input:
        csvs
    output:
        1 csv
    threads:
        16
    shell:
        cat

