#!/usr/bin/env python
# -*- coding: utf-8
#
# For usage, type: python average_images.py -h
# TODO: include this script inside processing script instead
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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)

spinal_levels = [2,3,4,5,6,7,8]
reg_type = ['PAM50_rootlets','Rootlets_reg', 'Discs_reg']
REG_COLOR = {'PAM50_rootlets':'green','Rootlets_reg':'blue', 'Discs_reg':'orange'}
COLOR_LEVEL = {2:'#8EFF00',3:'#00FF43',4:'#00FFF8',5:'#004FFF',6:'#8200FF',7:'#FF00D2',8:'#FF0009'}

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


def create_boxplot(df_overlap, output_path):
    # create boxplot with seaborn (x-ax is spinal level, y-ax is overlap):
    plt.figure()
    sns.set_style("darkgrid")

    # Create the boxplot with hue based on sex
    ax = sns.boxplot(x='level', y='overlap', data=df_overlap, hue='sex', boxprops=dict(alpha=.5), legend=False)
    colors = [["#4374B3"], ["darkorange"]]
    idx = 0
    for sex in df_overlap['sex'].unique():
        if sex == 'M':
            df_overlap['level'] = df_overlap['level'] - 0.2
        else:
            df_overlap['level'] = df_overlap['level'] + 0.4

        sns.set_palette(sns.color_palette(colors[idx]))
        sns.scatterplot(x='spinal_level', y='overlap', size='age', data=df_overlap[df_overlap['sex'] == sex],
                        sizes=(10, 150), hue='sex', edgecolor='w', legend=False)
        sns.set_style("darkgrid")
        idx += 1

    ax.set_title('Overlap between rootlets and PAM50 rootlets by method')
    ax.set_ylabel('Overlap [%]')
    ax.set_xlabel('Spinal level')
    ax.set_xticklabels(df_overlap['level'].unique())
    ax.set_yticks(range(0, 101, 10))

    lines = [
        Line2D([0], [0], color='#4374b3', linewidth=3, linestyle='-', label='Male'),
        Line2D([0], [0], color='#ff8c00', linewidth=3, linestyle='-', label='Female')
    ]
    age_circles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=5, label='Lower Age'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='Higher Age')
    ]
    handles = lines + age_circles

    # Create the legend
    plt.legend(handles=handles, loc='lower left', ncol=2, fontsize=10, handlelength=1)

    plt.savefig(f'{output_path}/boxplot_overlap_rootlets_vertebrae.svg', dpi=300)
    plt.show()


def create_coverage_plot(df, df_discs, df_pam50, path_out):
    fig = plt.figure(figsize=(5, 8))
    #ax = sns.boxplot(x='type_reg', y='start', data=df, hue='level', boxprops=dict(alpha=.5), legend=True)
    ax = fig.add_subplot()
        # Get the distance from PMJ and height of spinal level
    #for reg_type in TODO: add reg type
    i = 0

    for d in [df_pam50, df, df_discs]:
        for level in spinal_levels:
            reg = reg_type[i]
            start = d.loc[df["level"]==level,'start'].mean()
            start_std = d.loc[df["level"]==level,'start'].std()
            height = d.loc[df["level"]==level,'len'].mean()
            stop_std = d.loc[df["level"]==level,'stop'].std()
            print(level, start,'+/-', start_std, height,'+/-', stop_std)
            ax.add_patch(
                patches.Rectangle(
                    (0.9+i/3, start),      # (x,y)
                    0.2,            # width
                    height,         # height
                    facecolor=COLOR_LEVEL[level],#REG_COLOR[reg],
                    alpha=0.5,
                    edgecolor='none'
                ))
            ax.add_patch(
                patches.Rectangle(
                    (0.9+i/3, start-start_std),      # (x,y)
                    0.2,            # width
                    height+start_std+stop_std,         # height
                    facecolor=COLOR_LEVEL[level],#REG_COLOR[reg],
                    alpha=0.2,
                    edgecolor='none'
                ))

            # Add Level number ot each rectangle
            ax.text(
                0.9+i/3+0.007,     # x
                start,       # y
                'C'+str(level),
                horizontalalignment='left',
                verticalalignment='bottom',
                fontsize=8,
                color='white',
                path_effects=[pe.withStroke(linewidth=1, foreground='black')]
            )
            # Add mid slide number
            ax.text(
                0.9+i/3+0.21,     # x
                start+height/2,       # y
                int(start+height/2),
                horizontalalignment='left',
                verticalalignment='center',
                fontsize=8,
                color='black'#,
                #path_effects=[pe.withStroke(linewidth=1, foreground='black')]
            )
                            # Add mean value
            ax.plot(
                [0.9+i/3, 0.9+i/3+ 0.2],
                [start+height/2, start+height/2],
                color='black',
                linewidth=1,
                alpha=0.5,
                linestyle='dashed'
            )
        i+=1
    ax.set_xlim(0.8, 1.9)
    ax.set_ylim(min(df['stop'].min(), df['start'].min())*0.99,
                    max(df['stop'].max(), df['start'].max())*1.01)
    #ax.set_xticks(reg_type)
    ax.set_xticks([1,1/3+1,2/3+1])
    ax.set_xticklabels(reg_type)
    ax.set_ylabel('Slice (I-->S)', fontsize=14)
    ax.set_yticks(range(730, 990, 30))
    ax.grid(axis='y', alpha=0.2)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(f'{path_out}/boxplot_overlap_rootlets_vertebrae.png', dpi=300)


def get_start_stop_spinal_levels(file, type_reg=None, d_PAM50=None):
    logger.info(f"File: {file}")
    file_nib = nib.load(file)
    file_data = np.array(file_nib.get_fdata())
    subject = file.split("/")[-4]
    dict_list = []
    for spinal_level in spinal_levels:
        z_data_level_index = np.unique(np.where(file_data == spinal_level)[-1])
        z_data_level_index.sort()
        if len(z_data_level_index) > 0:
            start = z_data_level_index[0]
            stop = z_data_level_index[-1]
            length = len(z_data_level_index)
        else:
            start = np.nan
            stop = np.nan
            length = np.nan

        d = {
            'subject': subject,
            'level': spinal_level,
            'start': start,
            'stop': stop,
            'len': length,
            'type_reg': type_reg
        }
        #print("Levels:", spinal_level, "start: ", start, "stop: ", stop, "length: ", length)
        if d_PAM50 is not None:
            # Calculate the overlap directly using max and min functions:
            start_pam50 = d_PAM50.loc[d_PAM50['level']==spinal_level, ['start']].values[0][0]
            end_pam50 = d_PAM50.loc[d_PAM50['level']==spinal_level, ['stop']].values[0][0]

            overlap_start = max(start, start_pam50)
            overlap_end = min(stop, end_pam50)

            # Check if there is an overlap:
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                total_rootlets_length = stop - start
                overlap = (overlap_length / total_rootlets_length) * 100
            else:
                overlap = 0
            #print("Overlap", overlap, " %")
            d['overlap'] = overlap
        dict_list.append(d)
    return dict_list



def main():
    parser = get_parser()
    args = parser.parse_args()
    path_data = args.path_data
    path_file = args.path_file
    path_out = os.path.abspath(args.path_out)
    # Create output directory if does not exist
    if not os.path.exists(path_out):
        os.makedirs(path_out)
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
    dir_list = [f for f in glob.glob(path_data+'/*/anat/' + path_file + '/*rootlets_2template.nii.gz')]
    dir_list_discs = [f for f in glob.glob(path_data+'/*/anat/' + "reg_discs" + '/*rootlets_2template.nii.gz')]
    print("Number of files:, ", len(dir_list))
    # Filter out entries containing any substring from subjects_to_remove
    dir_list_excluded = [entry for entry in dir_list if not any(subject in entry for subject in dict_exclude_subj)]
    dir_list_discs_excluded = [entry for entry in dir_list_discs if not any(subject in entry for subject in dict_exclude_subj)]
    print("Number of files after excluding:, ", len(dir_list_excluded))
    #input_files = " ".join(dir_list_excluded[1:])  # put in strings for the command

    # Add PAM50 template
    path_sct = os.environ.get('SCT_DIR')
    path_rootlets_PAM50 = path_sct +'/data/PAM50/template/PAM50_rootlets.nii.gz'
    df_PAM50 = pd.DataFrame(get_start_stop_spinal_levels(path_rootlets_PAM50))
    print(df_PAM50)
    
    dict_list = []
    # Loop through rootlets file for reg_rootlets
    for file in dir_list_excluded:
        d = get_start_stop_spinal_levels(file, type_reg="reg_rootelts", d_PAM50=df_PAM50)
        dict_list.extend(d)
    #Loop through rootlets file for reg_discs
    dict_list_dict = []
    for file in dir_list_discs_excluded:
         d = get_start_stop_spinal_levels(file, type_reg="reg_discs", d_PAM50=df_PAM50)
         dict_list_dict.extend(d)

    df_coverage = pd.DataFrame(dict_list)
    df_coverage_discs = pd.DataFrame(dict_list_dict)
    # TODO: save in csv file
    df_coverage.to_csv(os.path.join(path_out,'coverage_rootlets_reg.csv'), index=False)
    df_coverage_discs.to_csv(os.path.join(path_out,'coverage_discs_reg.csv'), index=False)
    # Create boxplot:
    create_coverage_plot(df_coverage, df_coverage_discs, df_PAM50, path_out)

if __name__ == '__main__':
    main()
