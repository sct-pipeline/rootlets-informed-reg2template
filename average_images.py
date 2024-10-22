#!/usr/bin/env python
# -*- coding: utf-8
# TODO
#
# For usage, type: python average_images.py -h

# Authors: Sandrine Bédard


import argparse
import numpy as np
import os
import logging
import glob
import sys
import yaml
import subprocess
import nibabel as nib

# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)


def get_parser():
    parser = argparse.ArgumentParser(
        description="Averages loops through folder to find images and averages them together.")
    parser.add_argument('-path-data', required=True, type=str,
                        help="Path of data_processed folder.")
    parser.add_argument('-path-file', required=False, type=str,
                        default='reg_rootlets',
                        help="Output warping field filename.")
    parser.add_argument('-path-out', required=False, type=str,
                        help="Output warping field filename.")
    parser.add_argument('-exclude', required=False, type=str,
                        help="Exclude list. Ex: exclude.yml")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    path_data = args.path_data
    path_file = args.path_file
    path_out = args.path_out
    if args.exclude is not None: 
        # Check if input yml file exists
        if os.path.isfile(args.exclude):
            fname_yml = args.exclude
        else:
            sys.exit("ERROR: Input yml file {} does not exist or path is wrong.".format(args.exclude))
        with open(fname_yml, 'r') as stream:
            try:
                dict_exclude_subj = yaml.safe_load(stream)
                logger.info(f"Excluded subjects: {dict_exclude_subj}")
            except yaml.YAMLError as exc:
                logger.error(exc)
    else:
        # Initialize empty dict if no config yml file is passed
        dict_exclude_subj = dict()
    print(path_data+'/*/anat/' + path_file + '/anat2template.nii.gz')
    dir_list = [f for f in glob.glob(path_data+'/*/anat/' + path_file + '/anat2template.nii.gz')]
    print("Number of files:, ", len(dir_list))
    # Filter out entries containing any substring from subjects_to_remove
    dir_list_excluded = [entry for entry in dir_list if not any(subject in entry for subject in dict_exclude_subj)]
    print("Number of files after excluding:, ", len(dir_list_excluded))
    #input_files = " ".join(dir_list_excluded[1:])  # put in strings for the command
    file_nib = nib.load(dir_list_excluded[0])
    file_data = np.array(file_nib.get_fdata())
    sum_data = np.zeros(shape=file_data.shape)
    for file in dir_list_excluded:
        file_nib = nib.load(file)
        sum_data +=np.array(file_nib.get_fdata())
        print(file)
    mean_data = sum_data / len(dir_list_excluded)
    nii_mean= nib.Nifti1Image(mean_data, file_nib.affine)
    fname_out_levels = 'mean_' + path_file + '.nii.gz'
    nib.save(nii_mean, os.path.join(path_out, fname_out_levels))

    #cmd = f'sct_maths -i {dir_list_excluded[0]} -add {input_files} -o {os.path.join(path_out,'sum_' + path_file + '.nii.gz')}'
    #subprocess.run(cmd, shell=True)
    #cmd = f'sct_maths -i {os.path.join(path_out,'sum_' + path_file + '.nii.gz')} -div {str(len(dir_list_excluded))} -o {os.path.join(path_out,'mean_' + path_file + '.nii.gz')}'
    #subprocess.run(cmd, shell=True)
    
if __name__ == '__main__':
    main()
