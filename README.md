# Rootlets-based registration to PAM50 template
Spinal rootlets-informed registration to the template
## Installation
### Dependencies
* Spinal cord Toolbox [v7.1](https://github.com/spinalcordtoolbox/spinalcordtoolbox/releases/tag/7.1)
* Python 3.9
### Datasets
* Spine Generic Multi-subject dataset version [r20250310](https://github.com/spine-generic/data-multi-subject/releases/tag/r20250310)
* Spine Generic Single-subject dataset version [r20211115](https://github.com/spine-generic/data-single-subject/releases/tag/r20211115)
## Installation
Clone GitHub repository
```
git clone git@github.com:sct-pipeline/rootlets-informed-reg2template.git
cd rootlets-informed-reg2template/processing_script
```

Using SCT virtual environment
```
source ${SCT_DIR}/python/etc/profile.d/conda.sh
conda activate venv_sct
```
## Inter-subject Validation
Run batch **rootlets levels** and disc-based spinal cord cross-sectional area.
```
  cd $SCT_DIR
  git checkout sb/5111-csa-ref-rootlets
  cd $home/rootlets-informed-reg2template/processing_script
  sct_run_batch -jobs -1 -path-data <PATH_DATA> -script process_spine_generic_rootlets_csa.sh -path-out <PATH_OUT>
```
Run batch **spinal reg levels** and disc-based spinal cord cross-sectional area.
```
  cd $SCT_DIR
  git checkout sb/5111-csa-ref-rootlets-fromPAM50
  cd $home/rootlets-informed-reg2template/processing_script
  sct_run_batch -jobs -1 -path-data <PATH_DATA> -script process_spine_generic_spinallevel_csa.sh -path-out <PATH_OUT>
```

## Intra-subject Validation
Run batch **spinal reg levels** and disc-based spinal cord cross-sectional area on spine generic single subject dataset. It also computes the neck angle between C1-C2 and C5-C6 intervertebral discs.
```
  sct_run_batch -jobs -1 -path-data <PATH_DATA> -script /processing_script/process_spine_generic_single_spinal_csa.sh -path-out <PATH_OUT>
```

## Generating figures
Run the command once for Multi-subject, once for single subject.

```
python compare_methods.py compare_methods.py [-h] -path-SC-1 PATH_SC_1 -path-SC-2 PATH_SC_2 [-name-SC-1 NAME_SC_1]
                          [-name-SC-2 NAME_SC_2]  [-vertlevel VERTLEVEL]
                          [-ref-subject REF_SUBJECT] [-exclude EXCLUDE] [-angle ANGLE] [-path-out PATH_OUT]
Example command:

python compare_methods.py -path-SC-1 /home/ge.polymtl.ca/sabeda/processed_data/spine_generic_csa_rootlets_interseg_2026-06-18/results/rootlets_PAM50 -name-SC-1 rootlets -name-SC-2 spinal_level -path-SC-2 /home/ge.polymtl.ca/sabeda/GRAMES/process_results/spine_generic_csa_rootlets_2026-01-20-discs+spinallevel/results/spinal_level_PAM50  -exclude exclude_spine_generic_csa_analysis.yml -path-out /home/ge.polymtl.ca/sabeda/processed_data/spine_generic_csa_rootlets_interseg_2026-06-18/results/figures_spinal_level_rootlets_2026-06-19 -participant-file $SCT_DIR/data/PAM50_normalized_metrics/participants.tsv
```

