Data analysis of Oxford Nanopore Direct RNA seq  sequence specific experiment on HeLa cells

## Snakemake Rules:
### Signal for Individual Reads
**seq2mv:** Plots raw signal around specified position on reference transcript for one given read; uses _seq2mv_direct_RNA.py_

**Squigualiser:** Executes [squigualiser](https://github.com/hiruna72/squigualiser) to show raw signal to given read (html output)

### Signal Features for all Reads

**signal_summary:** Extracts signal features of all reads aligning to given reference positions (1337, 1842, 429) and saves them in one big .npz file, uses _signal_summary.py_

**Boxplot_Signal:** _boxplot_signal.py_ takes .npz file with generated features and makes boxplots for different given positions and features



### Other Rules:
**slice_bam:** creates bamfiles containing only a single read --> use for testing purposes (boxplot for one read should be similar to single read plot + squigualizer plot)
