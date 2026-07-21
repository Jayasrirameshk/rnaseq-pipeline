# importing libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.decomposition import PCA


# reading files
res_lum_vs_bas = pd.read_csv("results/deseq2_lum_vs_bas.csv", index_col = 0 )
res_preg_vs_vir = pd.read_csv("results/deseq2_preg_vs_vir.csv", index_col = 0 )
res_lac_vs_vir = pd.read_csv("results/deseq2_lac_vs_vir.csv", index_col = 0 )
res_lac_vs_preg = pd.read_csv("results/deseq2_lac_vs_preg.csv", index_col = 0 )
vst_count = pd.read_csv("results/vst_transformation.csv", index_col = 0)

# assigning significance
res_lum_vs_bas["significant"] = np.where((res_lum_vs_bas["padj"] < 0.05) & (abs(res_lum_vs_bas["log2FoldChange"]) > 1), "Significant", "Not Significant")
res_preg_vs_vir["significant"] = np.where((res_preg_vs_vir["padj"] < 0.05) & (abs(res_preg_vs_vir["log2FoldChange"]) > 1), "Significant", "Not Significant")
res_lac_vs_vir["significant"] = np.where((res_lac_vs_vir["padj"] < 0.05) & (abs(res_lac_vs_vir["log2FoldChange"]) > 1), "Significant", "Not Significant")
res_lac_vs_preg["significant"] = np.where((res_lac_vs_preg["padj"] < 0.05) & (abs(res_lac_vs_preg["log2FoldChange"]) > 1), "Significant", "Not Significant")

# volcano plots
plt.figure()
sns.scatterplot(data = res_lum_vs_bas, x = "log2FoldChange", y = -np.log10(res_lum_vs_bas['pvalue']), hue = "significant", palette = {"Significant": "red", "Not Significant" :"grey"})
plt.title("Volcano Plot for Luminal vs Basal")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(p-value)")
plt.savefig("results/figures/volcano_lum_vs_bas.png")
plt.close()

plt.figure()
sns.scatterplot(data = res_preg_vs_vir, x = "log2FoldChange", y = -np.log10(res_preg_vs_vir['pvalue']), hue = "significant", palette = {"Significant": "red", "Not Significant" :"grey"})
plt.title("Volcano Plot for Pregnant vs Virgin")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(p-value)")
plt.savefig("results/figures/volcano_preg_vs_vir.png")
plt.close()

plt.figure()
sns.scatterplot(data = res_lac_vs_vir, x = "log2FoldChange", y = -np.log10(res_lac_vs_vir['pvalue']), hue = "significant", palette = {"Significant": "red", "Not Significant" :"grey"})
plt.title("Volcano Plot for Lactate vs Virgin")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(p-value)")
plt.savefig("results/figures/volcano_lac_vs_vir.png")
plt.close()

plt.figure()
sns.scatterplot(data = res_lac_vs_preg, x = "log2FoldChange", y = -np.log10(res_lac_vs_preg['pvalue']), hue = "significant", palette = {"Significant": "red", "Not Significant" :"grey"})
plt.title("Volcano Plot for Lactate vs Pregnant")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(p-value)")
plt.savefig("results/figures/volcano_lac_vs_preg.png")
plt.close()

# prepping data for pca plots
vst_trans = vst_count.T
pca = PCA(n_components = 2)
pca_res = pca.fit_transform(vst_trans)
pc1_var = round(pca.explained_variance_ratio_[0] * 100, 1)
pc2_var = round(pca.explained_variance_ratio_[1] * 100, 1)
pca_df = pd.DataFrame(pca_res, columns=["PC1", "PC2"], index=vst_trans.index)
metadata = pd.read_csv("results/metadata.csv", index_col = 0)
pca_df["celltype"] = metadata["celltype"]
pca_df["condition"] = metadata["condition"]

# pca plots
plt.figure()
sns.scatterplot(data = pca_df, x = "PC1", y = "PC2", hue = "celltype", style = "condition" )
plt.title("PCA Plot")
plt.xlabel(f"PC1 ({pc1_var}%)")
plt.ylabel(f"PC2 ({pc2_var}%)")
plt.savefig("results/figures/pca_plot.png")
plt.close()

# prepping data for heatmap
var_gen = vst_count.var(axis=1)
top_50 = var_gen.nlargest(50).index
data_hm = vst_count.loc[top_50]
data_hm.columns = [col[:7] for col in data_hm.columns]


# heatmaps
sns.clustermap(data_hm, cmap = "RdBu_r", standard_scale=(10,12))
plt.savefig("results/figures/heatmap.png")
plt.close()