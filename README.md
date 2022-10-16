# HiTEA
The High Throughput Experimentation Analyzer (HiTEA) modules as described in King-Smith et al.

1. The workflow begins with the cleanup of a dataset as laid out in Example_Dataset_Cleanup.ipynb
2. Tanimoto Heatmaps are generated as laid out in Tanimoto_Heatmap.ipynb
3. Minor cleanup and one-hot encoding is performed as described in the Buchwald.R / Ullmann.R / Hetero_hydrogenation.R / Homo_hydrogenation.R
Note that these files are a compilation of numerous modules. The user should manually run the modules of interest. There are commented headings within the R files to assist in understanding which set of code and functions a user may be interested in running.
3. Random Forest variable importances are calculated as laid out in Random_Forest_Importances.ipynb
4. ANOVA-Tukey is then performed as described in Buchwald.R / Ullmann.R / Hetero_hydrogenation.R / Homo_hydrogenation.R
5. PCA of the ligands is performed as described in Example_Ligand_PCA.ipynb

# Important Note on Datasets
The datasets will not be released until the manuscript is accepted.
