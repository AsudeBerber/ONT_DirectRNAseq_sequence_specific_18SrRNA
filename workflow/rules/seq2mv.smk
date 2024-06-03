rule seq2mv:
    input: 
     "thisfile.doesntexist"
    output:
     ""
    conda:
        "../envs/seq2mv.yaml"
    shell:
     ""
     
