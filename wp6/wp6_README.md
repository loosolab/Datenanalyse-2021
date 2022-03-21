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

## Testdata
For Testdata ask @selinaLa or @oKoch
