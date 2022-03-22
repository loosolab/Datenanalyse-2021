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

###  Build tf_comb documentation:
If you want to have a look into the documentation of tfcomb, you have to do the following steps in your tfcomb repository.
1. ```make html```
2. Under **_builds/** you find the documentation. Open index.html with browser. 

## How to navigate:
You can find all notebooks/code of wp6 in the folder 'Datenanalyse-2021/wp6':
- ./Datenanalyse-2021/wp6/**analyse**: Here you will find our jupyter-notebooks with the tfcomb analysis and the imports of wp1/2/3/5 - Data
- ./Datenanalyse-2021/wp6/**testdata**: Here you find the used genome-files and jaspar-file for tfcomb, testdata for some examples of tfcomb and the other wp´s 

### ./Datenanalyse-2021/wp6/**analyse**:
  - **wp1_tf_cooccurences.ipynb**: Interface notebook for importing and analysing wp1´s - Data, details are in the notebook/ wiki
  - **wp2_tf_cooccurences.ipynb**: Interface notebook for importing and analysing wp2´s - Data, details are in the notebook/ wiki
  - **wp3_tf_co_occurences.ipynb**: Interface notebook for importing and analysing wp3´s - Data, details are in the notebook/ wiki
  - **wp5_new_motif_co_occurences.ipynb**: Interface notebook, merge new tf-motif´s of wp5 to an existing jaspar-motif file (run wp1/2/3 interface notbooks again with the adapted jaspar-file to get market-basket analysis with the new motifs).
  - **sl_enhancer_analyses_....ipynb**: contains example analysis to learn tfcomb for enhancer testdata
  - **sl_enhancer_analysis_combination_orientation_distance_....ipynb**: contains the combination of orientation and distance analysis of tfcomb for data of the wp´s 1/2/3
  - **sl_final_analysis_differenceindistanceofsameTFpairs_....ipynb**: contains analysis for differences in the binding distance of tfpairs in a tissue.
  - **sl_final_analysis_top50orientation-preferredbindingdistance....ipynb**:contains analysis for preffered binding distances of tf-pairs in a tissue.
  - **./results/**: contains the results of our analysis (e.g. tfcomb´s - market basket analysis, tfcomb´s - differential analysis, interesting/special tf-co-occurences)
     - ./results/**workpackage(e.g.wp1,2,3)/main/\<tissue>/\<cluster>.pkl**: contains tf-comb object´s with a executed market basket analysis for each cluster in the investigated tissue. They can be read-in again as a tfcomb-CombObj for further investigation.
     - ./results/**workpackage(e.g.wp1,2,3)/diff_analysis\<tissue>/\<analysis>.pkl**: contains a tf-comb diffObj with a differential analysis that is done for all cluster in a tissue, to find differences between the cluster of a tissue.
     - ./results/**workpackage(e.g.wp1,2,3)/answers**: contains Dataframes or csv with interesting tf-co-occurences, that are answers to the biological questions of wp6 
     
