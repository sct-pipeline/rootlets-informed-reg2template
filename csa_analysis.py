#!/usr/bin/env python
# -*- coding: utf-8

# Analyses CSA morphometry of using rootlets vs disc-based registration
# Example command: python csa_analysis.py -i ~/processed_data/results_spine_generic_csa/ -o ~/processed_data/results_spine_generic_csa/results_spine_generic_csa_2024-11-21 -exclude ~/code/rootlets-informed-reg2template/processing_script/exclude_spine_generic.yml
# Author: Sandrine Bédard


import os
import logging
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import scipy.stats as stats
import yaml
from scipy.signal import find_peaks


FNAME_LOG = 'log_stats.txt'

# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)

METRICS = ['MEAN(area)', 'MEAN(diameter_AP)', 'MEAN(diameter_RL)', 'MEAN(eccentricity)',
           'MEAN(solidity)']


METRICS_TO_YLIM = {
    'MEAN(diameter_AP)': (4, 9.3),
    'MEAN(area)': (30, 95),
    'MEAN(diameter_RL)': (8.5, 16),
    'MEAN(eccentricity)': (0.6, 0.95),
    'MEAN(solidity)': (0.912, 0.999),
    'std_smoothed_normalized_area': (0, 0.2),
    'smoothed_normalized_area': (0.9, 1.15)
}


METRIC_TO_AXIS = {
    'MEAN(diameter_AP)': 'AP Diameter [mm]',
    'MEAN(area)': 'Cross-Sectional Area [mm²]',
    'MEAN(diameter_RL)': 'Transverse Diameter [mm]',
    'MEAN(eccentricity)': 'Eccentricity [a.u.]',
    'MEAN(solidity)': 'Solidity [%]',
    'std_smoothed_normalized_area': 'STD (normalized CSA)',
    'smoothed_normalized_area': 'Normalized CSA'
}


PALETTE = {
    'sex': {'M': 'blue', 'F': 'red'},
    'group': {'rootlet': '#e31a1c', 'disc': 'blue'}
    }


LABELS_FONT_SIZE = 14
TICKS_FONT_SIZE = 12


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i-folder",
                        required=True,
                        type=str,
                        help="Results folder of spinal cord preprocessing")
    parser.add_argument("-o-folder",
                        type=str,
                        required=True,
                        help="Folder to write results")
    parser.add_argument('-ref-subject', required=False, type=str, default='sub-amu01',
                        help="Name of a reference subject ot use to get discs levels")

    parser.add_argument("-exclude",
                        type=str,
                        required=False,
                        help="exclude list")
    return parser


def format_pvalue(p_value, alpha=0.001, decimal_places=3, include_space=True, include_equal=True):
    """
    Format p-value.
    If the p-value is lower than alpha, format it to "<0.001", otherwise, round it to three decimals

    :param p_value: input p-value as a float
    :param alpha: significance level
    :param decimal_places: number of decimal places the p-value will be rounded
    :param include_space: include space or not (e.g., ' = 0.06')
    :param include_equal: include equal sign ('=') to the p-value (e.g., '=0.06') or not (e.g., '0.06')
    :return: p_value: the formatted p-value (e.g., '<0.05') as a str
    """
    if include_space:
        space = ' '
    else:
        space = ''

    # If the p-value is lower than alpha, return '<alpha' (e.g., <0.001)
    if p_value < alpha:
        p_value = space + "<" + space + str(alpha)
    # If the p-value is greater than alpha, round it number of decimals specified by decimal_places
    else:
        if include_equal:
            p_value = space + '=' + space + str(round(p_value, decimal_places))
        else:
            p_value = space + str(round(p_value, decimal_places))

    return p_value


def get_vert_indices(df, vertlevel='VertLevel'):
    """
    Get indices of slices corresponding to mid-vertebrae
    Args:
        df (pd.dataFrame): dataframe with CSA values
    Returns:
        vert (pd.Series): vertebrae levels across slices
        ind_vert (np.array): indices of slices corresponding to the beginning of each level (=intervertebral disc)
        ind_vert_mid (np.array): indices of slices corresponding to mid-levels
    """
    # Get vert levels for one certain subject
    vert = df[(df['participant_id'] == ref) & (df['group'] == 'disc')][vertlevel] # 'sub-amu01' TODO: add argument for example subject
    # Get indexes of where array changes value
    ind_vert = vert.diff()[vert.diff() != 0].index.values
    # Get the beginning of C1
    ind_vert = np.append(ind_vert, vert.index.values[-1])
    ind_vert_mid = []
    # Get indexes of mid-vertebrae
    for i in range(len(ind_vert)-1):
        ind_vert_mid.append(int(ind_vert[i:i+2].mean()))

    return vert, ind_vert, ind_vert_mid


def compare_metrics_across_group(df, perlevel=False, metric_chosen=None):
    """
    Compute Wilcoxon rank-sum tests between males and females for each metric.
    """

    logger.info("")

    for metric in METRICS:
        logger.info(f"\n{metric}")
        if metric_chosen:
            metric = metric_chosen
        if perlevel:
            slices_HC = df[df['group'] == 'rootlet'].groupby(['VertLevel'])[metric].mean()
            slices_HC_STD = df[df['group'] == 'rootlet'].groupby(['VertLevel'])[metric].std()
            slices_CR = df[df['group'] == 'disc'].groupby(['VertLevel'])[metric].mean()
            slices_CR_STD = df[df['group'] == 'disc'].groupby(['VertLevel'])[metric].std()
            logger.info(f'Mean {metric} for rootlet: {slices_HC}')
            logger.info(f'STD {metric} for rootlet: {slices_HC_STD}')
            logger.info(f'Mean {metric} for disc: {slices_CR}')
            logger.info(f'STD {metric} for disc: {slices_CR_STD}')
        else:

            # Get mean values for each slice
            slices_HC = df[df['group'] == 'rootlet'].groupby(['Slice (I->S)'])[metric].mean()
            slices_CR = df[df['group'] == 'disc'].groupby(['Slice (I->S)'])[metric].mean()

        # Run normality test
        stat, pval = stats.shapiro(slices_HC)
        logger.info(f'Normality test HC: p-value{format_pvalue(pval)}')
        stat, pval = stats.shapiro(slices_CR)
        logger.info(f'Normality test CR: p-value{format_pvalue(pval)}')
        # Run Wilcoxon rank-sum test (groups are independent)
        #from statsmodels.sandbox.stats.multicomp import multipletests
        #stat, pval = stats.ranksums(x=slices_HC, y=slices_CR)
        #p_adjusted = multipletests(pval, method='bonferroni')
        #print(p_adjusted)
        logger.info(f'{metric}: Wilcoxon rank-sum test between HC and CR: p-value{format_pvalue(pval)}')
        if metric_chosen:
            break


def read_t2w_pam50(file, group, exclude_list=None):
    df = pd.read_csv(file)
    df['participant_id'] = (df['Filename'].str.split('/').str[-4]).str.replace('_T2w_seg.nii.gz', '')
    df['group'] = group
    df = df[['participant_id', 'group', 'VertLevel', 'Slice (I->S)', 'MEAN(area)', 'MEAN(diameter_AP)', 'MEAN(diameter_RL)', 'MEAN(eccentricity)',
             'MEAN(solidity)']].drop(0)
    logger.info('Number of subjects BEFORE exclusion:')
    logger.info(len(np.unique(df['participant_id'])))
    if exclude_list:
        for subject in exclude_list:
            logger.info(f'dropping {subject}')
            df = df.drop(df[df['participant_id'] == subject].index, axis=0)
    logger.info('Number of subjects AFTER exclusion:')
    logger.info(len(np.unique(df['participant_id'])))
    return df


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def plot_ind_sub(df, group, metric, path_out, filename, hue='participant_id', y_min=0.8, y_max=1.2, peak_df=None):
    fig_size = (7, 6)
    font_size = LABELS_FONT_SIZE

    # Plot individual subject lines
    plt.figure()
    fig, ax = plt.subplots(figsize=fig_size)
    sns.lineplot(ax=ax, data=df.loc[df['group'] == group], x="Slice (I->S)", y=metric, hue=hue, linewidth=1, zorder=1, alpha=0.8)
    if peak_df is not None:
        sns.scatterplot(data=peak_df.loc[peak_df['group'] == group], x='Slice (I->S)', hue=hue, y=metric, style='participant_id', s=30, ax=ax, markers='o', zorder=2, edgecolor='black')
    ax.set_ylim(y_min, y_max)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(775, xmax - 15)
    ymin, ymax = ax.get_ylim()
    ax.get_legend().remove()

    # Get indices of slices corresponding vertebral levels
    vert, ind_vert, ind_vert_mid = get_vert_indices(df, vertlevel='VertLevel')
    for idx, x in enumerate(ind_vert[1:-1]):
        ax.axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)
    for idx, x in enumerate(ind_vert_mid, 0):
        if vert[x] > 7:
            level = 'T' + str(vert[x] - 7)
        else:
            level = 'C' + str(int(vert[x]))
        ax.text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                verticalalignment='bottom', color='black', fontsize=font_size)

    ax.invert_xaxis()
    ax.set_axisbelow(True)
    ax.tick_params(axis='both', which='major', labelsize=TICKS_FONT_SIZE)
    ax.set_ylabel("Normalized CSA", fontsize=font_size)
    ax.set_xlabel('Axial Slice #', fontsize=font_size)

    path_filename = os.path.join(path_out, filename)
    plt.savefig(path_filename, dpi=300, bbox_inches='tight')
    logger.info('Figure saved: ' + path_filename)

    # Create density plot with peaks
    if peak_df is not None:
        plt.figure()
        fig_kde, ax = plt.subplots(figsize=fig_size)
        sns.kdeplot(data=peak_df, x='Slice (I->S)', y=metric, hue='group', fill=True, common_norm=True, hue_order=['disc','rootlet'], alpha=0.7, ax=ax, palette=PALETTE['group'])
        ax.get_legend().remove()
        ax.set_xlabel('Axial Slice #', fontsize=font_size)
        ax.set_ylabel('Density', fontsize=font_size)
        ax.tick_params(axis='both', which='major', labelsize=TICKS_FONT_SIZE)
        ax.set_xlim(775, xmax - 15)
        #ymin, ymax = ax.get_ylim()
        ax.set_ylim(y_min, y_max)
        for idx, x in enumerate(ind_vert[1:-1]):
            ax.axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)
        for idx, x in enumerate(ind_vert_mid, 0):
            if vert[x] < 7:
                level = 'C' + str(int(vert[x]))
                ax.text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                        verticalalignment='bottom', color='black', fontsize=font_size)

        ax.invert_xaxis()
        plt.tight_layout()
        plt.savefig(os.path.join(path_out, "density.png"), dpi=300, bbox_inches='tight')
        logger.info(f'Figure saved: {os.path.join(path_out, "density.png")}')




def create_lineplot(df, metric=METRICS, hue=None, filename=None):
    """
    Create lineplot for individual metrics per vertebral levels.
    Note: we are ploting slices not levels to avoid averaging across levels.
    Args:
        df (pd.dataFrame): dataframe with metric values
        hue (str): column name of the dataframe to use for grouping; if None, no grouping is applied
    """

    #mpl.rcParams['font.family'] = 'Arial'

    fig, axes = plt.subplots(1, 5, figsize=(25, 4))
    axs = axes.ravel()

    # Loop across metrics
    for index, metric in enumerate(metric):
        # Note: we are ploting slices not levels to avoid averaging across levels
        if hue == 'sex' or hue == 'group':
            sns.lineplot(ax=axs[index], x="Slice (I->S)", y=metric, data=df, errorbar='sd', hue=hue, linewidth=2,
                         palette=PALETTE[hue])
            if index == 0:
                axs[index].legend(loc='upper right', fontsize=TICKS_FONT_SIZE)
            else:
                axs[index].get_legend().remove()
        else:
            sns.lineplot(ax=axs[index], x="Slice (I->S)", y=metric, data=df, errorbar='sd', hue=hue, linewidth=2)

        axs[index].set_ylim(METRICS_TO_YLIM[metric][0], METRICS_TO_YLIM[metric][1])
        ymin, ymax = axs[index].get_ylim()

        # Add labels
        axs[index].set_ylabel(METRIC_TO_AXIS[metric], fontsize=LABELS_FONT_SIZE)
        axs[index].set_xlabel('Axial Slice #', fontsize=LABELS_FONT_SIZE)
        # Increase xticks and yticks font size
        axs[index].tick_params(axis='both', which='major', labelsize=TICKS_FONT_SIZE)

        # Remove spines
        axs[index].spines['right'].set_visible(False)
        axs[index].spines['left'].set_visible(False)
        axs[index].spines['top'].set_visible(False)
        axs[index].spines['bottom'].set_visible(True)

        # Get indices of slices corresponding vertebral levels
        vert, ind_vert, ind_vert_mid = get_vert_indices(df)
        # Insert a vertical line for each intervertebral disc
        for idx, x in enumerate(ind_vert[:-1]):
            axs[index].axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)
        # Insert a text label for each vertebral level
        for idx, x in enumerate(ind_vert_mid, 0):
            # Deal with T1 label (C8 -> T1)
            if vert[x] > 7:
                level = 'T' + str(int(vert[x]) - 7)
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)
            else:
                level = 'C' + str(int(vert[x]))
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)

        # Invert x-axis
        axs[index].invert_xaxis()
        # Add only horizontal grid lines
        axs[index].yaxis.grid(True)
        # Move grid to background (i.e. behind other elements)
        axs[index].set_axisbelow(True)

    # Save figure
    if filename is None:
        if hue:
            filename = 'lineplot_per' + hue + '.png'
        else:
            filename = 'lineplot.png'
    plt.savefig(filename, dpi=500, bbox_inches='tight')
    logger.info('Figure saved: ' + filename)


def detect_peaks(data, y_col, x_col, height=None, prominence=None, distance=None, width=None):
    """
    Detect peaks in a specified column of a pandas DataFrame.
    
    Parameters:
        data (pd.DataFrame): Input DataFrame sorted by x-axis (e.g., Slice).
        y_col (str): Column name to detect peaks in (e.g., smoothed values).
        x_col (str): Column corresponding to x-axis (e.g., Slice).
        height (float): Minimum height of peaks.
        prominence (float): Minimum prominence of peaks.
        distance (int): Minimum distance between peaks (in terms of data points).

    Returns:
        pd.DataFrame: A DataFrame containing detected peaks with x and y values.
    """
    peaks_results = []

    # Group by 'group' and 'participant_id' to detect peaks separately for each
    for (group, participant), group_df in data.groupby(['group', 'participant_id']):
        y = group_df[y_col].values
        x = group_df[x_col].values

        # Detect peaks
        peaks, properties = find_peaks(y, height=height, prominence=prominence, distance=distance, width=width)

        if len(peaks) > 0:
            # If multiple peaks, choose the most prominent peak
            max_prominence_idx = np.argmax(properties['prominences'])
            peak_idx = peaks[max_prominence_idx]
        else:
            # If no peaks are detected, fall back to the maximum value
            peak_idx = np.argmax(y)

        # Store results
        peaks_results.append({
            'group': group,
            'participant_id': participant,
            x_col: x[peak_idx],
            y_col: y[peak_idx]
        })

    # Convert results to a DataFrame
    peaks_df = pd.DataFrame(peaks_results)
    return peaks_df


def normalize_csa(df):
    target_levels = [2.0, 3.0]  # VertLevel values to target
    slice_range = 5  # Number of slices above and below
    # Filter rows corresponding to the target vertebral levels
    level_filtered = df[df['VertLevel'].isin(target_levels)]

    # Step 2: Identify the slice ranges for each participant and vert level
    # Merge the original DataFrame with the level_filtered DataFrame to map participant_id and Slice (I->S) ranges
    merged = df.merge(
        level_filtered[['participant_id', 'Slice (I->S)']],
        on='participant_id',
        suffixes=('', '_target')
    )

    # Step 3: Filter rows within ±slice_range of the target Slice (I->S)
    slices_of_interest = merged[
        (merged['Slice (I->S)'] >= merged['Slice (I->S)_target'] - slice_range) &
        (merged['Slice (I->S)'] <= merged['Slice (I->S)_target'] + slice_range)
    ]
    # Step 4: Compute the mean of MEAN(area) for each participant
    result = (
        slices_of_interest
        .groupby('participant_id')['MEAN(area)']
        .mean()
        .reset_index()
    )

    # Rename columns for clarity
    result.columns = ['participant_id', 'average_mean_area']

    # Merge the average CSA back into the original dataframe
    df_with_avg = df.merge(result, on='participant_id', how='left')
    # Normalize MEAN(area) by dividing by the average CSA
    df_with_avg['normalized_mean_area'] = df_with_avg['MEAN(area)'] / df_with_avg['average_mean_area']

    df_with_avg = df_with_avg.sort_values(by=['participant_id', 'Slice (I->S)'])

    # Step 3: Apply the smoothing function to the normalized_mean_area column
    box_pts = 55  # Smoothing window size
    df_with_avg['smoothed_normalized_area'] = (
                df_with_avg.groupby('participant_id')['normalized_mean_area']
                .transform(lambda x: smooth(x, box_pts))
                )
    # Compute STD of normalized csa
    df_with_avg['std_smoothed_normalized_area'] = df_with_avg.groupby('Slice (I->S)')['smoothed_normalized_area'].std()
    return df_with_avg


def main():

    args = get_parser().parse_args()
    # Get input argments
    input_folder = os.path.abspath(args.i_folder)
    global ref
    ref = args.ref_subject
    output_folder = args.o_folder
    # Create output folder if does not exist.
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    os.chdir(output_folder)
    # Dump log file there
    if os.path.exists(FNAME_LOG):
        os.remove(FNAME_LOG)
    fh = logging.FileHandler(os.path.join(FNAME_LOG))
    logging.root.addHandler(fh)

    # Create a list with subjects to exclude if input .yml config file is passed
    if args.exclude is not None:
        # Check if input yml file exists
        if os.path.isfile(args.exclude):
            fname_yml = args.exclude
        else:
            sys.exit("ERROR: Input yml file {} does not exist or path is wrong.".format(args.exclude))
        with open(fname_yml, 'r') as stream:
            try:
                exclude = list(yaml.safe_load(stream))
            except yaml.YAMLError as exc:
                logger.error(exc)
    else:
        exclude = []
    logger.info(exclude)

    # Analyse T2w perslice
    #################################################################
    logger.info('\nAnalysing T2w CSA perslice in PAM50 anatomical dimension')
    filename_rootlets = os.path.join(input_folder, "csa_rootlets_PAM50.csv")
    filename_discs = os.path.join(input_folder, "csa_discs_PAM50.csv")

    df_rootlets = read_t2w_pam50(filename_rootlets, 'rootlet', exclude_list=exclude).dropna(axis=0)
    df_discs = read_t2w_pam50(filename_discs, 'disc', exclude_list=exclude).dropna(axis=0)
    df_rootlets = normalize_csa(df_rootlets)
    df_discs = normalize_csa(df_discs)
    df_all = pd.concat([df_rootlets, df_discs], ignore_index=True)
    df_all = df_all[df_all['VertLevel'] < 7]
    df_all = df_all[df_all['VertLevel'] > 1]

    logger.info("Detect peak of cervical enlargement")
    peaks_df = detect_peaks(
               df_all,
               y_col='smoothed_normalized_area',
               x_col='Slice (I->S)',
               height=None,       # Set minimum height (optional)
               prominence=0.05,    # Set minimum prominence to filter significant peaks
               distance=None,         # Minimum distance between peaks
               width=None
    )
    print(peaks_df['Slice (I->S)'].max())
    # Print the 10 maximum values and related subjects
    top_10_peaks = peaks_df.nlargest(10, 'Slice (I->S)')
    print("Top 10 maximum values and related subjects:")
    print(top_10_peaks[['participant_id', 'smoothed_normalized_area', 'Slice (I->S)']])
    peaks_df.to_csv(os.path.join(output_folder, 'peaks.csv'))
    logger.info('Number of participants:')
    logger.info(len(np.unique(df_all['participant_id'])))
    create_lineplot(df_all, hue='group')
    plot_ind_sub(df_all, group='rootlet', metric='smoothed_normalized_area', path_out=output_folder, filename='csa_persubject_normalized_rootlet.png', peak_df=peaks_df)
    plot_ind_sub(df_all, group='disc', metric='smoothed_normalized_area', path_out=output_folder, filename='csa_persubject_normalized_disc.png', peak_df=peaks_df)

    plot_ind_sub(df_all, group='rootlet', metric='MEAN(area)', path_out=output_folder, filename='csa_persubject_rootlet.png', y_min=40, y_max=110)
    plot_ind_sub(df_all, group='disc', metric='MEAN(area)', path_out=output_folder, filename='csa_persubject_disc.png', y_min=40, y_max=110)
    create_lineplot(df_all, metric=['std_smoothed_normalized_area'], hue='group', filename='lineplot_std_normalized_csa.png')
    create_lineplot(df_all, metric=['smoothed_normalized_area'], hue='group', filename='lineplot_normalized_csa.png')
    df_all = df_all[df_all['VertLevel'] < 7]
    df_all = df_all[df_all['VertLevel'] > 1]

    mean_csa_rootlets = df_all[df_all['group'] == 'rootlet']['MEAN(area)'].mean()
    std_csa_rootlets = df_all[df_all['group'] == 'rootlet']['MEAN(area)'].std()

    mean_csa_rootlets_norm = df_all[df_all['group'] == 'rootlet']['normalized_mean_area'].mean()
    std_csa_rootlets_norm = df_all[df_all['group'] == 'rootlet']['normalized_mean_area'].std()

    cov_rootlets = (std_csa_rootlets/mean_csa_rootlets)*100
    logger.info('COV\n')
    logger.info(cov_rootlets)
    logger.info(f'Mean CSA for rootlets reg between C2 and C7: {mean_csa_rootlets} +- {std_csa_rootlets}')
    logger.info(f'Mean norm CSA for rootlets reg between C2 and C7: {mean_csa_rootlets_norm} +- {std_csa_rootlets_norm}')

    # Compute metrics for discs-based reg
    mean_csa_discs = df_all[df_all['group'] == 'disc']['MEAN(area)'].mean()
    mean_csa_discs_norm = df_all[df_all['group'] == 'disc']['normalized_mean_area'].mean()

    std_csa_discs = df_all[df_all['group'] == 'disc']['MEAN(area)'].std()
    std_csa_discs_norm = df_all[df_all['group'] == 'disc']['normalized_mean_area'].std()

    cov_discs = (std_csa_discs/mean_csa_discs)*100
    # normalized_mean_area
    logger.info('COV\n')
    logger.info(cov_discs)

    logger.info(f'Mean CSA for disc reg between C2 and C7: {mean_csa_discs} +- {std_csa_discs}')
    logger.info(f'Mean CSA for disc reg between C2 and C7: {mean_csa_discs_norm} +- {std_csa_discs_norm}')


if __name__ == "__main__":
    main()
