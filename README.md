# HiTEA
The High Throughput Experimentation Analyzer (HiTEA) modules as described in King-Smith et al.

## Workflow
1. The workflow begins with the cleanup of a dataset as laid out in Example_Dataset_Cleanup.ipynb
2. Tanimoto Heatmaps are generated as laid out in Tanimoto_Heatmap.ipynb
3. Minor cleanup and one-hot encoding is performed as described in the Buchwald.R / Ullmann.R / Hetero_hydrogenation.R / Homo_hydrogenation.R
Note that these files are a compilation of numerous modules. The user should manually run the modules of interest. There are commented headings within the R files to assist in understanding which set of code and functions a user may be interested in running.
3. Random Forest variable importances are calculated as laid out in Random_Forest_Importances.ipynb
4. ANOVA-Tukey is then performed as described in Buchwald.R / Ullmann.R / Hetero_hydrogenation.R / Homo_hydrogenation.R
5. PCA of the ligands is performed as described in Example_Ligand_PCA.ipynb

## Important Note on Datasets
The datasets will not be released until the manuscript is accepted.

## Dependencies
All code is run on python 3.7 and R 3.4.4.
##### Python Dependencies:
* scikitlearn==1.0.1
* numpy==1.21.1
* pandas==1.3.4
* rdkit==2020.09.01
* matplotlib==3.3.4
* mols2grid==0.1.0
* [adjustText](https://github.com/Phlya/adjustText)
* IPython==7.29.0
* re==2.2.1
##### R Dependencies:
* data.table==1.12.2
* mltools==0.3.5
* stringi==1.4.3
