rule merge_pod5:
    input:
        "resources/pod5/{dir}/"
    output:
        "resources/pod5/{dir}/pod5_merge/pod5.pod5"
    threads:
        16
    shell:
        "pod5 merge {input}*.pod5 -o {output}"


rule sort_pod5:
    input: 
        "resources/pod5/{dir}/pod5_concat/.pod5"
    output:
        "resources/pod5/{dir}/pod5_concat/sorted.pod5"
    shell:
        "sort"
