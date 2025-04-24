# Rootlets-based registration to PAM50 template
Spinal rootlets-informed registration to the template
## Installation
### Dependencies
* Spinal cord Toolbox [v6.5](https://github.com/spinalcordtoolbox/spinalcordtoolbox/releases/tag/6.5)
* Python 3.9
### Datasets
* Spine Generic Multi-subject dataset version [r20250310](https://github.com/spine-generic/data-multi-subject/releases/tag/r20250310)
* Spinal Cord Head Positions dataset version [v1.1.1](https://openneuro.org/datasets/ds004507/versions/1.1.1)
## Instalation
Clone Github repository
```
git clone git@github.com:sct-pipeline/rootlets-informed-reg2template.git
cd rootlets-informed-reg2template
```

Using SCT virtual environement
```
source ${SCT_DIR}/python/etc/profile.d/conda.sh
conda activate venv_sct
```
## Inter-subject Validation
Run batch regsitration to template using rootlets and disc-based template regsitration. This will also compute spianl cord cross-sectional area in template space
```
  sct_run_batch -jobs -1 -path-data <PATH_DATA> -script /processing_script/process_spine_generic.sh -path-out <PATH_OUT>
```
Create average image of rootlets-based registration:
```
python average_images.py -path-data <PATH_DATA_PROCESSED> -path-file reg_rootlets -path-out <PATH_RESULTS> -exclude processing_script/exclude_spine_generic.yml 
```
For disc-based registration
```
python average_images.py -path-data <PATH_DATA_PROCESSED> -path-file reg_discs -path-out <PATH_RESULTS> -exclude processing_script/exclude_spine_generic.yml 
```

### Morphometric analysis

```
python csa_analysis.py -i <PATH_RESULTS> -o <PATH_RESULTS>/results_spine_generic_csa_2024-11-21 -exclude processing_script/exclude_spine_generic.yml
```
## Single-subject Validation

```
  sct_run_batch -jobs -1 -path-data <PATH_DATA> -script /processing_script/process_pmj_data.sh -path-out <PATH_OUT>
```

## Analyse Functional Data

```
python analyse_func_ridge.py -i1 <rootlets_smooth_3.1_FE/thresh_zstat1_med.nii.gz> -i2 <discs_smooth_3.1_FE/thresh_zstat1_med.nii.gz> -o <PATH_OUT>
```
