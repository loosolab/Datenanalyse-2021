# WP6
WP6 is searching for transcription factor (TF)-co-occurences in the clustered single cell ATAC-Seq-Data with the python library tfcomb (https://github.com/loosolab/TF-COMB). Tfcomb is using a market basket analysis for finding (TF)-co-occurences. For further details of tfcomb and the Project background, please have look at our WIKI (https://github.com/loosolab/Datenanalyse-2021/wiki/WP6). 
The ReadMe handles the setup for WP6 and how to use the jupyter-notebooks for investigating WP1/WP2/WP3´s data.  

## Setup for WP6´s conda environment:
We use jupyter-notebooks for our analysis. To use our notebooks you have to do the following steps for setting up your environment.
Tfcomb is still under development, because of this you have to install it from source.

1. Install miniconda (https://docs.conda.io/en/latest/miniconda.html#installing)
2. Checkout tf_comb (https://github.com/loosolab/TF-COMB)
3. cd tfcomb
4. Create conda environment: ```conda create -n tfcomb_env --file required_packages.txt```
5. Activate Conda environment: ```conda activate tfcomb_env ```
6. Install jupyter-notebooks: ```conda install ipykernel```
7. Register kernel for jupyter notebook: ```python -m ipykernel install --name tfcomb_env --user```
8. In the tf_comb repository do: ```pip install .```  , this installs tfcomb into your envrionment.

##  Build tf_comb documentation:
If you want to have a look into the documentation of tfcomb, you have to do the following steps in your tfcomb repository.
1. ```make html```
2. Under **_builds/** you find the documentation. Open index.html with browser. 

## How to navigate:
You can find all notebooks/code of wp6 in the folder 'Datenanalyse-2021/wp6':
- ./Datenanalyse-2021/wp6/**analyse**: Here you will find our jupyter-notebooks with the tfcomb analysis and the imports of wp1/2/3/5 - Data
- ./Datenanalyse-2021/wp6/**testdata**: Here you find the used genome-files and jaspar-file for tfcomb, testdata for some examples of tfcomb and the other wp´s 

### ./Datenanalyse-2021/wp6/**analyse**:
  - **./results/**: contains the results of our analysis (e.g. tfcomb´s - market basket analysis, tfcomb´s - differential analysis, interesting/special tf-co-occurences)
     - ./results/**workpackage(e.g.wp1)/main/\<tissue>/\<cluster>.pkl**: contains tf-comb object´s with a exectued market basket analysis for each cluster in the inverstigated tissue. They can be read-in again as a tfcomb-CombObj for further investigation.
     - ./results/**workpackage(e.g.wp1)/diff_analysis\<tissue>/\<analysis>.pkl**: contains a tf-comb diffObj with a differential analysis that is done for all cluster in a tissue, to find differences between the cluster of a tissue.
     - ./results/**workpackage(e.g.wp1)/answers**: contains Dataframes or csv with interesting tf-co-occurences, that are answers to the biological questions of wp6
