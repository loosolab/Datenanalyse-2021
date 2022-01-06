# WP6
WP6 is searching for co-occurences of transcription factors (TF) with a market basket analyses. For this we use the tfcomb library.

## Example Setup conda environment:
1. Checkout tf_comb (https://github.com/loosolab/TF-COMB)
2. cd tfcomb
3. Create conda environment: ```conda create -n tfcomb_env --file required_packages.txt```
4. Activate Conda environment: ```conda activate tfcomb_env ```
5. ```conda install ipykernel```
6. register kernel for jupyter notebook: ```python -m ipykernel install --name tfcomb_env --user```
7. In the tf_comb repository do: ```pip install .```

##  Build tf_comb documentation:
1. ```make html```
2. Under **_builds/** you find the documentation. Open index.html with browser. 

## Testdata
For Testdata ask @selinaLa or @oKoch