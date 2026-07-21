# Loading relevant libraries
library(DESeq2)
library(yaml)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)

# Reading the config file
config = read_yaml("config/config.yaml")

# Unpacking the config into variables
counts_file = config$counts
results_dir = config$results_dir
samples_info = config$samples

# Setting up figure directory
figures_dir = file.path(results_dir,"figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Reading the raw data file 
raw_counts = read.delim(counts_file, row.names = "EntrezGeneID", check.names = FALSE)

# Preparation of data for DESeq2 analysis  
counts = raw_counts[,colnames(raw_counts) != "Length"]

sample_names = names(samples_info)
celltype = sapply(samples_info, function(x) x$celltype)
condition = sapply(samples_info, function(x)x$condition)

col_Data = data.frame (celltype = factor(celltype, levels = c("basal","luminal")), condition = factor (condition, levels = c("virgin","pregnant","lactate")),row.names = sample_names)
counts = counts [,sample_names]

# Creation of the DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = counts, colData = col_Data, design = ~ condition + celltype)

# Filtering low count genes
dds = dds[rowSums(counts(dds)) >= 10,]

# Running DESeq
dds = DESeq(dds)

# Extracting results
res_celltype = results(dds, contrast = c("celltype", "luminal", "basal"))
res_preg_vs_virgin = results(dds, contrast = c("condition", "pregnant", "virgin"))
res_lact_vs_virgin = results(dds, contrast = c("condition", "lactate", "virgin"))
res_lact_vs_preg = results(dds, contrast = c("condition", "lactate", "pregnant"))

# Saving results as a csv file
res_celltype_df = as.data.frame(res_celltype)
res_lact_vs_preg_df = as.data.frame(res_lact_vs_preg)
res_lact_vs_virgin_df = as.data.frame(res_lact_vs_virgin)
res_preg_vs_virgin_df = as.data.frame(res_preg_vs_virgin)

write.csv(res_celltype_df, file = file.path(results_dir, "deseq2_lum_vs_bas.csv"), row.names = TRUE)
write.csv(res_lact_vs_preg_df, file = file.path(results_dir, "deseq2_lac_vs_preg.csv"), row.names = TRUE)
write.csv(res_lact_vs_virgin_df, file = file.path(results_dir, "deseq2_lac_vs_vir.csv"), row.names = TRUE)
write.csv(res_preg_vs_virgin_df, file = file.path(results_dir, "deseq2_preg_vs_vir.csv"), row.names = TRUE)

write.csv(col_Data, file = file.path(results_dir, "metadata.csv"), row.names = TRUE)

# Saving output files for visualization
saveRDS(dds, file = file.path(results_dir, "deseq2_res.rds"))

# VST transformation
vsd = vst(dds, blind = TRUE)
write.csv(assay(vsd), file = file.path(results_dir, "vst_transformation.csv"), row.names = TRUE)

# Adding gene annotation
gene_symbols = mapIds(org.Mm.eg.db,
                      keys = rownames(res_celltype_df),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")

res_celltype_df$gene_symbol = gene_symbols[rownames(res_celltype_df)]
res_preg_vs_virgin_df$gene_symbol = gene_symbols[rownames(res_preg_vs_virgin_df)]
res_lact_vs_virgin_df$gene_symbol = gene_symbols[rownames(res_lact_vs_virgin_df)]
res_lact_vs_preg_df$gene_symbol = gene_symbols[rownames(res_lact_vs_preg_df)]

write.csv(res_celltype_df, file = file.path(results_dir, "deseq2_lum_vs_bas.csv"), row.names = TRUE)
write.csv(res_preg_vs_virgin_df, file = file.path(results_dir, "deseq2_preg_vs_vir.csv"), row.names = TRUE)
write.csv(res_lact_vs_virgin_df, file = file.path(results_dir, "deseq2_lac_vs_vir.csv"), row.names = TRUE)
write.csv(res_lact_vs_preg_df, file = file.path(results_dir, "deseq2_lac_vs_preg.csv"), row.names = TRUE)