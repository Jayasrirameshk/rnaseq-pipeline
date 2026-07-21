configfile: "config/config.yaml"

rule all:
     input:
            "results/figures/volcano_lum_vs_bas.png",
            "results/figures/volcano_preg_vs_vir.png",
            "results/figures/volcano_lac_vs_vir.png",
            "results/figures/volcano_lac_vs_preg.png",
            "results/figures/pca_plot.png",
            "results/figures/heatmap.png"

rule deseq2:
     input:
          counts = config["counts"],
          config = "config/config.yaml"
     output:
           "results/deseq2_lum_vs_bas.csv",
           "results/deseq2_preg_vs_vir.csv",
           "results/deseq2_lac_vs_vir.csv",
           "results/deseq2_lac_vs_preg.csv",
           "results/vst_transformation.csv",
           "results/metadata.csv"
     script:
           "scripts/deseq2.R"

rule plots:
    input:
        "results/deseq2_lum_vs_bas.csv",
        "results/deseq2_preg_vs_vir.csv",
        "results/deseq2_lac_vs_vir.csv",
        "results/deseq2_lac_vs_preg.csv",
        "results/vst_transformation.csv",
        "results/metadata.csv"
    output:
        "results/figures/volcano_lum_vs_bas.png",
        "results/figures/volcano_preg_vs_vir.png",
        "results/figures/volcano_lac_vs_vir.png",
        "results/figures/volcano_lac_vs_preg.png",
        "results/figures/pca_plot.png",
        "results/figures/heatmap.png"
    script:
        "scripts/plotting.py"