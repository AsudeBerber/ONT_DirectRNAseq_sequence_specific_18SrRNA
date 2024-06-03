rule seq2mv_single_read:
    input: 
     "thisfile.doesntexist"
    output:
     "../../resources/signal/{bam_dir}/plots/{read_ids}/{read_ids}_{start}-{end}.png"
    params:
        read = "x"
        bam = 
        # change to pod5/p2i for p2i data
        pod5 = "pod5/p2s"
        base_pos = 
        bases_around = 
    conda:
        "../envs/seq2mv.yaml"
    shell:
     "python seq2mv_direct_RNA.py {wildcards.dir} RNAseq_test_noisy_correction_sorted {params.read} 300 350 --pod5_dir resources/{pod5}"
     
 