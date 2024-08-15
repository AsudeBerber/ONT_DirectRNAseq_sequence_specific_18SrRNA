include: "workflow/rules/dorado.smk"
include: "workflow/rules/get_transcriptome.smk"
include: "workflow/rules/minimap2_alignment.smk"
include: "workflow/rules/seq2mv.smk"
include: "workflow/rules/samtools.smk"
include: "workflow/rules/squigualiser.smk"
include: "workflow/rules/slice_bam.smk"
include: "workflow/rules/boxplot_signal.smk"
include: "workflow/rules/signal_summary.smk"

rule all:
    input:
        "resources/signal/p2s/plots/fc27240e-46d4-4f79-9aab-a61efb10860d/fc27240e-46d4-4f79-9aab-a61efb10860d_1842-pm8.svg",
        "resources/signal/p2s/plots/00009085-ae8c-45e3-aea9-b9c3f784433e/00009085-ae8c-45e3-aea9-b9c3f784433e_1842-pm8.svg",
        "resources/signal/p2s/signal_summary/CCG_window_21_p2s_aligned_sorted/1337_1842_430_event_1.svg",
        "resources/signal/squigualiser/fc27240e-46d4-4f79-9aab-a61efb10860d-PAW35875_9fd38647_68d05f77_357/p2s_aligned_1800-1850.html"
