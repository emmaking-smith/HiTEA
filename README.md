# HiTEA
The **Hi**gh **T**hroughput **E**xperimentation **A**nalyzer (HiTEA) modules as described in King-Smith *et al.*
Please note: Analysis for the manuscript was performed with the R modules (.R files in R_modules directory) however, python modules for ANOVA, Tukey, and Random Forest tests are also provided. Given the differences between implementations of 2-way ANOVA and Tukey HSD in R and python, some minor differences may be observed.

[![DOI](https://zenodo.org/badge/552294062.svg)](https://zenodo.org/badge/latestdoi/552294062)


## Installing Juypter Notebook
Some modules (.ipynb) files are jupyter notebook files and require installation of jupyter lab. See the [official website](https://jupyter.org/install) for details of installation and running a notebook on your specific system.

## Modules of HiTEA
HiTEA is made up of several modules that can be used on their own or in conjuction with one another.
Before any analysis can be performed, the dataset(s) must be cleaned. An example of how to accomplish this is laid out in `Example_Dataset_Cleanup.ipynb`.

### Variable Importances
Each of the four datasets described in the manuscript have their own specific module, denoted as `Random_Forest_XX.py` where `XX` is the type of reaction being analyzed. All random forest modules can be run with default parameters, however, specific arguments can be passed. For example:

`python Random_Forest_Buchwald.py -d {PATH_TO_YOUR_DATASET} -s {PATH_TO_DESIRED_SAVE_DIRECTORY} -t {SPECIFIC_YEAR_TO_ANALYZE} -n {REMOVE 0% YIELDING REACTIONS TRUE / FALSE}`

Specifiying a `-t` argument yields a temporal analysis of that specific year, and specifying `-n` as `True` yields the analysis of the dataset without its 0% yielding reactions. 

The output will be saved as a csv file in the specified directory. Note that for temporal analyses, the corresponding dataset is `XX_with_time.csv` where `XX` is the type of reaction needs to be used.

### ANOVA-Tukey
Similar to the variable importances, four specific modules for ANOVA-Tukey analysis are provided denoted as `ANOVA_Tukey_XX.py` where `XX` is the type of reaction being analyzed. Again, all modules can be run with default parameters, with temporal analysis or no 0% yielding reaction analysis performed with the relevant flags. For example:

`python ANOVA_Tukey_Buchwald.py -d {PATH_TO_YOUR_DATASET} -s {PATH_TO_DESIRED_SAVE_DIRECTORY} -t {SPECIFIC_YEAR_TO_ANALYZE} -n {REMOVE 0% YIELDING REACTIONS TRUE / FALSE}`

The output will be saved as a csv file in the specified directory. Note that for temporal analyses, the corresponding dataset is `XX_with_time.csv` where `XX` is the type of reaction needs to be used.

### Visualizing the Chemical Space
Two modules are provided for visualizing the chemical space. These are `Example_Ligand_PCA.ipynb` and `Tanimoto_Heatmap.ipynb` which visualize the ligand space and reactant space, respectively.

## Important Note on Datasets
The datasets will not be released until the manuscript is accepted. 

***Note For Reviewers:*** 
You should have received access to the full dataset as well as the cleaned Buchwald, Ullmann, Heterogeneous Hydrogenation, and Homogeneous Hydrogenation datasets. Please reach out to esk34@cam.ac.uk or alpha@alpha-lee.com if this is not the case.

## Dependencies
All code has been tested on python 3.7 and R 3.4.4.
#### Python Dependencies:
Install the following dependencies with: `pip install X` where `X` is the name as written below.
* scikit-learn
* numpy
* pandas
* matplotlib
* mols2grid
* adjustText (see [here](https://github.com/Phlya/adjustText) for GitHub repo)
* ipython
* statsmodels

Install RDKit with conda:

`conda install -c rdkit rdkit`

#### R Dependencies (if using R modules):
* data.table==1.12.2
* mltools==0.3.5
* stringi==1.4.3

## License:
MIT License
