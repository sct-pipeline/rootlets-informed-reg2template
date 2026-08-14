#
# Functions to plot morphometric metrics computed from normative database (spine-generic dataset in PAM50 space)
# perslice and vertebral levels
#
# Example usage:
#       python compare_methods.py
#       -path-HC $SCT_DIR/data/PAM50_normalized_metrics
#       -participant-file $SCT_DIR/data/PAM50_normalized_metrics/participants.tsv
# python compare_methods.py -path-SC-1 ~/code/PAM50-normalized-metrics -path-SC-2  ~/duke/temp/sebeda/spine_generic_canal_seg_2024-12-05/results/ -participant-file ~/code/PAM50-normalized-metrics/participants.tsv
# Author: Sandrine Bédard, Jan Valosek
#

import os
import argparse
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import yaml


FNAME_LOG = 'log_stats.txt'

# Initialize logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # default: logging.DEBUG, logging.INFO
hdlr = logging.StreamHandler(sys.stdout)
logging.root.addHandler(hdlr)

METRICS = ['MEAN(area)', 'MEAN(diameter_AP)', 'MEAN(diameter_RL)', 'MEAN(compression_ratio)', 'MEAN(eccentricity)',
           'MEAN(solidity)']

METRICS_DTYPE = {
    'MEAN(diameter_AP)': 'float64',
    'MEAN(area)': 'float64',
    'MEAN(diameter_RL)': 'float64',
    'MEAN(eccentricity)': 'float64',
    'MEAN(solidity)': 'float64',
    'aSCOR':'float64'
}

METRIC_TO_TITLE = {
    'MEAN(diameter_AP)': 'AP Diameter',
    'MEAN(area)': 'Cross-Sectional Area',
    'MEAN(diameter_RL)': 'Transverse Diameter',
    'MEAN(eccentricity)': 'Eccentricity',
    'MEAN(solidity)': 'Solidity',
    'MEAN(compression_ratio)': 'Compression Ratio',
    'aSCOR':'aSCOR',
}

METRIC_TO_AXIS = {
    'MEAN(diameter_AP)': 'AP Diameter [mm]',
    'MEAN(area)': 'Cross-Sectional Area [mm²]',
    'MEAN(diameter_RL)': 'Transverse Diameter [mm]',
    'MEAN(eccentricity)': 'Eccentricity [a.u.]',
    'MEAN(solidity)': 'Solidity [%]',
    'MEAN(compression_ratio)': 'Compression Ratio [a.u.]',
    'aSCOR':'aSCOR [%]',
}

DEMOGRAPHIC_TO_AXIS = {
    'age': 'Age [years]',
    'BMI': 'BMI [kg/m²]',
    'height': 'Height [cm]',
    'weight': 'Weight [kg]',
}

# ylim max offset (used for showing text)
METRICS_TO_YLIM_OFFSET = {
    'MEAN(diameter_AP)': 0.4,
    'MEAN(area)': 6,
    'MEAN(diameter_RL)': 0.7,
    'MEAN(eccentricity)': 0.03,
    'MEAN(solidity)': 1,
    'MEAN(compression_ratio)': 0.03,
    'aSCOR': 1,
}

# Set ylim to do not overlap horizontal grid with vertebrae labels
METRICS_TO_YLIM = {
    'MEAN(diameter_AP)': (5, 9.3), #(6,8.5), #TODO: use second value for canal
    'MEAN(area)': (30, 90), # (40, 85)
    'MEAN(diameter_RL)': (8, 14.5), #(15, 35),
    'MEAN(eccentricity)': (0.65, 0.9),
    'MEAN(solidity)': (91.2, 99.9),
    'MEAN(compression_ratio)': (0.41, 0.84),
    'aSCOR': (20, 50),
}

DISCS_DICT = {
    7: 'C7-T1',
    6: 'C6-C7',
    5: 'C5-C6',
    4: 'C4-C5',
    3: 'C3-C4',
    2: 'C2-C3',
    1: 'C1-C2'
}

MID_VERT_DICT = {
    8: 'T1',
    7: 'C7',
    6: 'C6',
    5: 'C5',
    4: 'C4',
    3: 'C3',
    2: 'C2',
    1: 'C1'
}

VENDORS = ['Siemens', 'Philips', 'GE']
AGE_DECADES = ['10-20', '21-30', '31-40', '41-50', '51-60']

LABELS_FONT_SIZE = 14
TICKS_FONT_SIZE = 12

COLORS_SEX = {
    'M': 'blue',
    'F': 'red'
    }

# To be same as spine-generic figures (https://github.com/spine-generic/spine-generic/blob/master/spinegeneric/cli/generate_figure.py#L114)
PALETTE = {
    'sex': {'M': 'blue', 'F': 'red'},
    'manufacturer': {'Siemens': 'green', 'Philips': 'dodgerblue', 'GE': 'black'},
    'age': {'10-20': 'blue', '21-30': 'green', '31-40': 'black', '41-50': 'red', '51-60': 'purple'},
    'method': {'DeepSegSC': 'blue', 'ContrastAgn2.5': 'red', 'rootlets': 'orange', 'discs': 'blue', 'spinal_level': 'red'}
    }
color_palette = {
        'sct_deepseg_sc': '#e78ac3',
        'contrast-agnostic_v2.0': '#66c2a5',
        'contrast-agnostic_v3.0': '#a6d854',
        'DWI':'#66c2a5',
        'MT-on':'#fc8d62',
        'GRE-T1w':'#8da0cb',
        'T1w':'#e78ac3',
        'T2*w':'#a6d854',
        'T2w': '#ffd92f',
        'rootlets': 'orange',
        'discs': 'blue',
        'spinal_level': 'red'
    }

def get_parser():
    parser = argparse.ArgumentParser(
        description="Plot morphometric metrics computed from normative database (spine-generic dataset in PAM50 "
                    "space) perslice and vertebral levels ")
    parser.add_argument('-path-SC-1', required=True, type=str,
                        help="Path to data of normative dataset computed perslice 1.")
    parser.add_argument('-path-SC-2', required=True, type=str,
                        help="Path to data of normative dataset computed perslice 2.")
    parser.add_argument('-name-SC-1', required=False, type=str, default='DeepSegSC',
                        help="Name of method/contrast or else of -path-SC-1")
    parser.add_argument('-name-SC-2', required=False, type=str, default='ContrastAgn2.5',
                        help="Name of method/contrast or else of -path-SC-2")
    parser.add_argument('-participant-file', required=False, type=str,
                        help="Path to participants.tsv file.")
    parser.add_argument('-vertlevel', required=False, type=str, default='1:9',
                        help="Vertebral levels to select (e.g., '2:8' for C2 to T1).")
    parser.add_argument('-ref-subject', required=False, type=str, default='sub-amu01',
                        help="Path to participants.tsv file.")
    parser.add_argument('-exclude', required=False, type=str,
                        help="Path to exclude.yml file.")
    parser.add_argument('-angle', required=False, type=str,
                        help="CSV file with neck angle.")
    parser.add_argument('-path-out', required=False, type=str, default='stats',
                        help="Output directory name.")
    # TODO: add exclude list option
    # TODO: add names of the methods
    return parser


def csv2dataFrame(filename):
    """
    Loads a .csv file and builds a pandas dataFrame of the data
    Args:
        filename (str): filename of the .csv file
    Returns:
        data (pd.dataFrame): pandas dataframe of the .csv file's data
    """
    data = pd.read_csv(filename)
    return data


def get_vert_indices(df):
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
    try:
        if ref in df[df['method'] == name_SC_1]['participant_id'].values:
            vert = df[(df['participant_id'] == ref) & (df['method'] == name_SC_1)]['VertLevel']
        elif ref in df[df['method'] == name_SC_2]['participant_id'].values:
            vert = df[(df['participant_id'] == ref) & (df['method'] == name_SC_2)]['VertLevel']
        else:
            raise KeyError(f"Reference subject {ref} not found in either method.")
    except KeyError:
       # logger.info('trying something else')
        vert = df[(df['participant_id'] == ref)]['VertLevel']
    # Get indexes of where array changes value
    ind_vert = vert.diff()[vert.diff() != 0].index.values
    # Get the beginning of C1
    ind_vert = np.append(ind_vert, vert.index.values[-1])
    ind_vert_mid = []
    # Get indexes of mid-vertebrae
    for i in range(len(ind_vert)-1):
        ind_vert_mid.append(int(ind_vert[i:i+2].mean()))

    return vert, ind_vert, ind_vert_mid


def create_lineplot(df, hue, path_out, show_cv=False, set_axis=True, filename=None):
    """
    Create lineplot for individual metrics per vertebral levels.
    Note: we are ploting slices not levels to avoid averaging across levels.
    Args:
        df (pd.dataFrame): dataframe with metric values
        hue (str): column name of the dataframe to use for grouping; if None, no grouping is applied
        path_out (str): path to output directory
        show_cv (bool): if True, include coefficient of variation for each vertebral level to the plot
    """

   # mpl.rcParams['font.family'] = 'Arial'

    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    axs = axes.ravel()

    # Loop across metrics
    for index, metric in enumerate(METRICS):
        # Note: we are ploting slices not levels to avoid averaging across levels
        if hue == 'sex' or hue == 'manufacturer' or hue == 'age':
            sns.lineplot(ax=axs[index], x="Slice (I->S)", y=metric, data=df, errorbar='sd', hue=hue, linewidth=2)#,
                         #palette=PALETTE[hue])
        else:
            sns.lineplot(ax=axs[index], x="Slice (I->S)", y=metric, data=df, errorbar='sd', hue=hue, palette=color_palette, linewidth=2)
            axs[index].set_xlim(722, 948)
            if index == 0:
                axs[index].legend(loc='lower center', bbox_to_anchor=(0.38, 0.1), fontsize=TICKS_FONT_SIZE)
            else:
                legend = axs[index].get_legend()
                if legend:
                    legend.remove()
        if set_axis:
            axs[index].set_ylim(METRICS_TO_YLIM[metric][0], METRICS_TO_YLIM[metric][1])
        else:
            axs[index].set_ylim(0.85, 1.2)
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
        for idx, x in enumerate(ind_vert[1:-1]):
            axs[index].axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)

        # Insert a text label for each vertebral level
        for idx, x in enumerate(ind_vert_mid, 0):
            if show_cv:
                cv = compute_cv(df[(df['VertLevel'] == vert[x])], metric)
            # Deal with T1 label (C8 -> T1)
            if vert[x] > 8:
                level = 'T' + str(vert[x] - 8)
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)
                # Show CV
                if show_cv:
                    axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymax-METRICS_TO_YLIM_OFFSET[metric],
                                    str(round(cv, 1)) + '%', horizontalalignment='center', verticalalignment='bottom',
                                    color='black')
            else:
                level = 'C' + str(vert[x])
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)
                # Show CV
                if show_cv:
                    axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymax-METRICS_TO_YLIM_OFFSET[metric],
                                    str(round(cv, 1)) + '%', horizontalalignment='center', verticalalignment='bottom',
                                    color='black')
            if show_cv:
                logger.info(f'{metric}, {level}, COV: {cv}')

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
    path_filename = os.path.join(path_out, filename)
    plt.savefig(path_filename, dpi=300, bbox_inches='tight')
    logger.info(f'Figure saved: {path_filename}')


def create_regplot(df, path_out, show_cv=False):
    """
    Plot data and a linear regression model fit. Slices in X and Coefficient of Variation (CoV) in Y.
    Args:
        df (pd.dataFrame): dataframe with metric values
        path_out (str): path to output directory
        show_cv (bool): if True, include coefficient of variation for each vertebral level to the plot
    """

    #mpl.rcParams['font.family'] = 'Arial'

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axs = axes.ravel()

    # Loop across metrics
    for index, metric in enumerate(METRICS):
        slices_list = []
        cv_list = []
        # Loop across slices
        for slice in df['Slice (I->S)'].unique():
            # Get metric value for each slice
            df_slice = df[df['Slice (I->S)'] == slice]
            cv_list.append(compute_cv(df_slice, metric))
            slices_list.append(slice)

        sns.regplot(ax=axs[index], x=slices_list, y=cv_list, scatter_kws={'alpha': 0.5})

        # Add labels
        axs[index].set_xlabel('Axial Slice #', fontsize=LABELS_FONT_SIZE)
        axs[index].set_ylabel(f'{METRIC_TO_TITLE[metric]} COV [%]', fontsize=LABELS_FONT_SIZE)
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
        for idx, x in enumerate(ind_vert[1:-1]):
            axs[index].axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)

        # Set the same y-axis limits across metrics
        axs[index].set_ylim([0, 18])

        # Place text box with mean COV value
        axs[index].text(.5, .94, 'COV: {}±{}%'.format(round(np.mean(cv_list), 1), round(np.std(cv_list), 1)),
                 horizontalalignment='center', verticalalignment='center', transform=axs[index].transAxes,
                 fontsize=TICKS_FONT_SIZE, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
        # Move the text box to the front
        axs[index].set_zorder(1)

        ymin, ymax = axs[index].get_ylim()
        # Insert a text label for each vertebral level
        for idx, x in enumerate(ind_vert_mid, 0):
            if show_cv:
                cv = compute_cv(df[(df['VertLevel'] == vert[x])], metric)
            # Deal with T1 label (C8 -> T1)
            if vert[x] > 7:
                level = 'T' + str(vert[x] - 7)
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)
                # Show CV
                if show_cv:
                    axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], 14.8,
                                    str(round(cv, 1)) + '%', horizontalalignment='center', verticalalignment='bottom',
                                    color='black')
            else:
                level = 'C' + str(vert[x])
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)
                # Show CV
                if show_cv:
                    axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], 14.8, str(round(cv, 1)) + '%',
                                    horizontalalignment='center', verticalalignment='bottom', color='black')
            if show_cv:
                logger.info(f'{metric}, {level}, COV: {cv}')

        # Invert x-axis
        axs[index].invert_xaxis()
        # Add only horizontal grid lines
        axs[index].yaxis.grid(True)
        # Move grid to background (i.e. behind other elements)
        axs[index].set_axisbelow(True)

    # Save figure
    filename = 'cov_scatterplot.png'
    path_filename = os.path.join(path_out, filename)
    plt.savefig(path_filename, dpi=300, bbox_inches='tight')
    logger.info(f'Figure saved: {path_filename}')


def create_regplot_per_methods(df, path_out):
    """
    Plot data and a linear regression model fit. Slices in X and Coefficient of Variation (CoV) in Y. Per sex.
    Args:
        df (pd.dataFrame): dataframe with metric values
        path_out (str): path to output directory
    """

    #mpl.rcParams['font.family'] = 'Arial'

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axs = axes.ravel()

    # Loop across metrics
    for index, metric in enumerate(METRICS):
        mean_cov = dict()
        std_cov = dict()
        # Loop across sex
        for sex in df['method'].unique():
            slices_list = []
            cv_list = []
            # Loop across slices
            for slice in df['Slice (I->S)'].unique():
                # Get metric value for each slice
                df_slice = df[(df['Slice (I->S)'] == slice) & (df['sex'] == sex)]
                cv_list.append(compute_cv(df_slice, metric))
                slices_list.append(slice)

            mean_cov[sex] = np.mean(cv_list)
            std_cov[sex] = np.std(cv_list)
            sns.regplot(ax=axs[index], x=slices_list, y=cv_list, label=sex, scatter_kws={'alpha': 0.5},
                        color=PALETTE['method'][sex])

        if index == 0:
            axs[index].legend(loc='upper right', fontsize=TICKS_FONT_SIZE)

        # Add labels
        axs[index].set_xlabel('Axial Slice #', fontsize=LABELS_FONT_SIZE)
        axs[index].set_ylabel(f'{METRIC_TO_TITLE[metric]} COV [%]', fontsize=LABELS_FONT_SIZE)
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
        for idx, x in enumerate(ind_vert[1:-1]):
            axs[index].axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)

        # Set the same y-axis limits across metrics
        axs[index].set_ylim([0, 18])

        # Place text box with COV values
        # Note: we invert xaxis, thus xmax is used for the left limit
        axs[index].text(.5, .90, 'F COV: {}±{}%\nM COV: {}±{}%'.format(round(mean_cov['F'], 1),
                                                                       round(std_cov['F'], 1),
                                                                       round(mean_cov['M'], 1),
                                                                       round(std_cov['M'], 1)),
                        horizontalalignment='center', verticalalignment='center', transform=axs[index].transAxes,
                        fontsize=TICKS_FONT_SIZE, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
        # Move the text box to the front
        axs[index].set_zorder(1)

        ymin, ymax = axs[index].get_ylim()
        # Insert a text label for each vertebral level
        for idx, x in enumerate(ind_vert_mid, 0):
            # Deal with T1 label (C8 -> T1)
            if vert[x] > 7:
                level = 'T' + str(vert[x] - 7)
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)
            else:
                level = 'C' + str(vert[x])
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)

        # Invert x-axis
        axs[index].invert_xaxis()
        # Add only horizontal grid lines
        axs[index].yaxis.grid(True)
        # Move grid to background (i.e. behind other elements)
        axs[index].set_axisbelow(True)

    # Save figure
    filename = 'cov_scatterplot_persex.png'
    path_filename = os.path.join(path_out, filename)
    plt.savefig(path_filename, dpi=300, bbox_inches='tight')
    logger.info(f'Figure saved: {path_filename}')


def compute_cv(df, metric):
    """
    Compute coefficient of variation (CV) of a given metric.
    Args:
        df (pd.dataFrame): dataframe with CSA values
        metric (str): column name of the dataframe to compute CV
    Returns:
        cv (float): coefficient of variation
    """
    cv = df[metric].std() / df[metric].mean()
    cv = cv * 100
    return cv


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


def compute_c2_c3_stats(df):
    """
    Compute mean and std from C2 and C3 levels across sex and compare females and males.
    """

    # Compute mean and std from C2 and C3 levels
    df_c2_c3 = df[(df['VertLevel'] == 2) | (df['VertLevel'] == 3)]
    c2_c3_persex = df_c2_c3.groupby('sex')['MEAN(area)'].agg([np.mean, np.std])
    logger.info(c2_c3_persex)

    # Compare C2-C3 CSA between females and males
    c2_c3_f = df_c2_c3[df_c2_c3['sex'] == 'F']['MEAN(area)']
    c2_c3_m = df_c2_c3[df_c2_c3['sex'] == 'M']['MEAN(area)']
    # Run normality test
    stat, pval = stats.normaltest(c2_c3_f)
    logger.info(f'Normality test C2-C3 females: p-value{format_pvalue(pval)}')
    stat, pval = stats.normaltest(c2_c3_m)
    logger.info(f'Normality test C2-C3 males: p-value{format_pvalue(pval)}')
    # Compute Mann-Whitney U test
    stat, pval = stats.mannwhitneyu(c2_c3_f, c2_c3_m)
    logger.info(f'Mann-Whitney U test between females and males: p-value{format_pvalue(pval)}')


def compare_metrics_across_methods(df):
    """
    Compute Wilcoxon rank-sum tests between methods for each metric.
    """

    logger.info("")

    for metric in METRICS:
        logger.info(f"\n{metric}")
        #restrict:    
        df_short = df[df['Slice (I->S)'] <= 948]
        df_short = df[df['Slice (I->S)'] >= 712]
        # Get mean values for each slice
        a = df_short[df_short['method'] == name_SC_1].groupby(['participant_id'])[metric].mean()
        b = df_short[df_short['method'] == name_SC_2].groupby(['participant_id'])[metric].mean()
        diff = a - b
        print("N =", len(diff))
        print("mean diff =", diff.mean())
        print("std diff =", diff.std(ddof=1))
        print(diff.describe())
        # Run normality test
        stat, pval = stats.shapiro(a-b)
        logger.info(f'Normality test: p-value{format_pvalue(pval)}')
        # Run Wilcoxon rank-sum test (groups are independent)
        
        if pval >= 0.05:
            stat, p = stats.ttest_rel(a, b)
            logger.info(f"{metric}: paired t-test: p={format_pvalue(p)}")
        else:
            stat, p = stats.wilcoxon(a, b)
            logger.info(f"{metric}: Wilcoxon signed-rank (paired): p={format_pvalue(p)}")
        
        # Repeat tests per vertebral level
        logger.info(f"  Per vertebral level:")
        for vert_level in sorted(df_short['VertLevel'].unique()):
            df_vert = df_short[df_short['VertLevel'] == vert_level]
            a_vert = df_vert[df_vert['method'] == name_SC_1].groupby(['participant_id'])[metric].mean()
            b_vert = df_vert[df_vert['method'] == name_SC_2].groupby(['participant_id'])[metric].mean()
            
            # Only test if both methods have data for this level
            common_ids_vert = a_vert.index.intersection(b_vert.index)
            if len(common_ids_vert) > 0:
                a_vert = a_vert.loc[common_ids_vert]
                b_vert = b_vert.loc[common_ids_vert]
                
                stat, pval_vert = stats.shapiro(a_vert - b_vert)
                if pval_vert >= 0.05:
                    stat, p_vert = stats.ttest_rel(a_vert, b_vert)
                    logger.info(f"    VertLevel {vert_level}: paired t-test: p={format_pvalue(p_vert)}")
                else:
                    stat, p_vert = stats.wilcoxon(a_vert, b_vert)
                    logger.info(f"    VertLevel {vert_level}: Wilcoxon signed-rank: p={format_pvalue(p_vert)}")
        
        mean_deepseg = df[df['method'] == name_SC_1][metric].mean()
        std_deepseg = df[df['method'] == name_SC_1][metric].std()
        logger.info(f'\n{name_SC_1}: mean +/- std = {mean_deepseg} +/- {std_deepseg}')

        mean_ca = df[df['method'] == name_SC_2][metric].mean()
        std_ca = df[df['method'] == name_SC_2][metric].std()
        logger.info(f'{name_SC_2}: mean +/- std = {mean_ca} +/- {std_ca}')
        
        # Compute average scaling factor:
        mean_scale = mean_ca/mean_deepseg
        logger.info(f'Scale: mean = {mean_scale}')



def compute_scaling_factor(df):
    df_scaling = pd.DataFrame()
    df_scaling['participant_id'] = df[df['method'] == name_SC_1]['participant_id']
    df_scaling['Slice (I->S)'] =  df[df['method'] == name_SC_1]['Slice (I->S)']
    df_scaling['VertLevel'] =  df[df['method'] == name_SC_1]['VertLevel']

    for metric in METRICS:
        logger.info(f"\n{metric}")
        slices_deepseg = df[df['method'] == name_SC_1].groupby(['Slice (I->S)', 'participant_id'])[metric].mean()
       # slices_deepseg = df[df['method'] == name_SC_1].groupby(['participant_id'])[metric].mean()
        slices_ca = df[df['method'] == name_SC_2].groupby(['Slice (I->S)', 'participant_id'])[metric].mean()
        scaling_factor = (slices_ca/slices_deepseg).unstack()
        scaling_factor = scaling_factor.reset_index().melt(value_name=metric, id_vars=['Slice (I->S)'], value_vars=np.unique(df_scaling['participant_id']))
        df_scaling = df_scaling.merge(scaling_factor, on=['Slice (I->S)', 'participant_id'])
        
        # Compute average scaling factor
        mean_scale = scaling_factor[metric].mean()
        std_scale = scaling_factor[metric].std()
        logger.info(f'Scaling factor: mean +/- std = {mean_scale} +/- {std_scale}')
    return df_scaling


def compute_normative_values(df, path_out):
    """
    Compute normative values for each metric and save them in a csv file
    Args:
        df:
        path_out:

    Returns:
    """
    for metric in METRICS:
        logger.info(f"\n{metric}")
        # Compute mean and std for each slice across all subjects
        slices_mean = df.groupby(['Slice (I->S)'])[metric].mean()
        slices_std = df.groupby(['Slice (I->S)'])[metric].std()

        # Get indices of slices corresponding vertebral levels
        vert, ind_vert, ind_vert_mid = get_vert_indices(df)

        d = []
        # Loop across intervertebral discs
        for x in reversed(ind_vert[1:-1]):
            slice_number = df.loc[x, 'Slice (I->S)']
            disc = DISCS_DICT[vert[x]]
            slice_mean = slices_mean.loc[slice_number]
            slice_std = slices_std.loc[slice_number]
            logger.info(f'Disc {disc}, slice {slice_number}: {round(slice_mean, 2)} ± {round(slice_std, 2)}')
            d.append(
                {
                    'Disc': disc,
                    'Slice': slice_number,
                    'Mean ± STD': f'{round(slice_mean, 2)} ± {round(slice_std, 2)}'
                }
            )

        fname_csv = os.path.join(path_out, metric + '_disc_normative_values.csv')
        pd.DataFrame(d).to_csv(fname_csv, index=False)
        logger.info(f'Created: {fname_csv}.\n')

        d = []
        # Loop across mid-vertebral slices
        for x in reversed(ind_vert_mid):
            slice_number = df.loc[x, 'Slice (I->S)']
            mid_level = MID_VERT_DICT[vert[x]]
            slice_mean = slices_mean.loc[slice_number]
            slice_std = slices_std.loc[slice_number]
            logger.info(f'Level {mid_level}, slice {slice_number}: {round(slice_mean, 2)} ± {round(slice_std, 2)}')
            d.append(
                {
                    'Level': mid_level,
                    'Slice': slice_number,
                    'Mean ± STD': f'{round(slice_mean, 2)} ± {round(slice_std, 2)}'
                }
            )

        fname_csv = os.path.join(path_out, metric + '_mid_level_normative_values.csv')
        pd.DataFrame(d).to_csv(fname_csv, index=False)
        logger.info(f'Created: {fname_csv}.\n')


def compute_normative_values_persex(df, path_out):
    """
    Compute normative values for each metric persex and save them in a csv file
    Args:
        df:
        path_out:

    Returns:
    """
    for metric in METRICS:
        logger.info(f"\n{metric}")
        # Compute mean and std for each slice across all subjects
        slices_mean = df.groupby(['Slice (I->S)', 'sex'])[metric].mean()
        slices_std = df.groupby(['Slice (I->S)', 'sex'])[metric].std()

        # Get indices of slices corresponding vertebral levels
        vert, ind_vert, ind_vert_mid = get_vert_indices(df)

        d = []
        # Loop across intervertebral discs
        for x in reversed(ind_vert[1:-1]):
            slice_number = df.loc[x, 'Slice (I->S)']
            disc = DISCS_DICT[vert[x]]
            # males
            slice_mean_M = slices_mean.loc[slice_number]['M']
            slice_std_M = slices_std.loc[slice_number]['M']
            # females
            slice_mean_F = slices_mean.loc[slice_number]['F']
            slice_std_F = slices_std.loc[slice_number]['F']
            logger.info(f'M, Disc {disc}, slice {slice_number}: {round(slice_mean_M, 2)} ± {round(slice_std_M, 2)}')
            logger.info(f'F, Disc {disc}, slice {slice_number}: {round(slice_mean_F, 2)} ± {round(slice_std_F, 2)}')
            d.append(
                {
                    'Disc': disc,
                    'Slice': slice_number,
                    'M, Mean ± STD': f'{round(slice_mean_M, 2)} ± {round(slice_std_M, 2)}',
                    'F, Mean ± STD': f'{round(slice_mean_F, 2)} ± {round(slice_std_F, 2)}'
                }
            )

        fname_csv = os.path.join(path_out, metric + '_disc_normative_values_persex.csv')

        pd.DataFrame(d).to_csv(fname_csv, index=False)
        logger.info(f'Created: {fname_csv}.\n')

        d = []
        # Loop across mid-vertebral slices
        for x in reversed(ind_vert_mid):
            slice_number = df.loc[x, 'Slice (I->S)']
            mid_level = MID_VERT_DICT[vert[x]]
            # males
            slice_mean_M = slices_mean.loc[slice_number]['M']
            slice_std_M = slices_std.loc[slice_number]['M']
            # females
            slice_mean_F = slices_mean.loc[slice_number]['F']
            slice_std_F = slices_std.loc[slice_number]['F']
            logger.info(f'M, Level {mid_level}, slice {slice_number}: {round(slice_mean_M, 2)} ± {round(slice_std_M, 2)}')
            logger.info(f'F, Level {mid_level}, slice {slice_number}: {round(slice_mean_F, 2)} ± {round(slice_std_F, 2)}')
            d.append(
                {
                    'Level': mid_level,
                    'Slice': slice_number,
                    'M, Mean ± STD': f'{round(slice_mean_M, 2)} ± {round(slice_std_M, 2)}',
                    'F, Mean ± STD': f'{round(slice_mean_F, 2)} ± {round(slice_std_F, 2)}'
                }
            )

        fname_csv = os.path.join(path_out, metric + '_mid_level_normative_values_persex.csv')
        pd.DataFrame(d).to_csv(fname_csv, index=False)
        logger.info(f'Created: {fname_csv}.\n')


def plot_metrics_relative_to_age(df, path_out_figures):
    """
    # Plot averaged metrics from C2-C3 levels as a function of age separated for sex
    Args:
        df:
        path_out_figures:

    Returns:

    """
    # Compute mean and std from C2 and C3 levels per sex
    df_c2_c3 = df[(df['VertLevel'] == 2) | (df['VertLevel'] == 3)]

    # Recode sex to 0 and 1 to all pd.mean() at the next step
    df_c2_c3['sex'].replace(['F', 'M'], [0, 1], inplace=True)

    for metric in METRICS:
        # Get mean values for each subject
        df_c2_c3_mean = df_c2_c3.groupby(['participant_id'])[[metric, 'age', 'sex']].mean()

        # Plot averaged metrics from C2-C3 levels as a function of age
        fig, ax = plt.subplots(figsize=(20, 6))
        sns.regplot(x='age', y=metric, data=df_c2_c3_mean[df_c2_c3_mean['sex'] == 1],
                    label='Male', color=COLORS_SEX['M'])
        sns.regplot(x='age', y=metric, data=df_c2_c3_mean[df_c2_c3_mean['sex'] == 0],
                    label='Female', color=COLORS_SEX['F'])
        ax.legend()
        ax.set_xlabel('Age (years)')
        ax.set_ylabel(METRIC_TO_AXIS[metric])
        ax.set_title(f'{METRIC_TO_TITLE[metric]} from C2-C3 levels as a function of age')

        # Save figure
        fname_fig = os.path.join(path_out_figures, metric + '_C2_C3_vs_age_persex.png')
        fig.savefig(fname_fig, bbox_inches='tight')
        logger.info(f'Created: {fname_fig}.\n')


def read_csv_files(path_HC, participant_file=None, exclude_list=None):
    # Initialize pandas dataframe where data across all subjects will be stored
    logger.info(f'Reading {path_HC}')
    df = pd.DataFrame()
    # Loop through .csv files of healthy controls
    for file in os.listdir(path_HC):
        if 'PAM50.csv' in file:
            # Read csv file as pandas dataframe for given subject
            df_subject = pd.read_csv(os.path.join(path_HC, file), dtype=METRICS_DTYPE)
            # Compute compression ratio (CR) as MEAN(diameter_AP) / MEAN(diameter_RL)
            df_subject['MEAN(compression_ratio)'] = df_subject['MEAN(diameter_AP)'] / df_subject['MEAN(diameter_RL)']

            # Concatenate DataFrame objects
            df = pd.concat([df, df_subject], axis=0, ignore_index=True)
    # Get sub-id (e.g., sub-amu01) from Filename column and insert it as a new column called participant_id
    # Subject ID is the first characters of the filename till slash
    df.insert(0, 'participant_id', df['Filename'].str.split('/').str[-3])
    # Get number of unique subjects (unique strings under Filename column)
    subjects = df['Filename'].unique() 
    # If a participants.tsv file is provided, insert columns sex, age and manufacturer from df_participants into df
    df_participants = pd.DataFrame()
    if participant_file:
        df_participants = pd.read_csv(participant_file, sep='\t')
        df = df.merge(df_participants[["age", "sex", "height", "weight", "manufacturer", "participant_id"]],
                    on='participant_id')
    # Exclude subjects if exclude_list is provided
    logger.info(f'Number of subjects before exclusion: {len(subjects)}')
    if exclude_list is not None:
        df = df[~df['participant_id'].isin(exclude_list)]
        subjects = df['Filename'].unique() 
        logger.info(f'Excluded subjects: {exclude_list}')
    # Print number of subjects
    logger.info(f'Number of subjects after exclusion: {str(len(subjects))}\n')
    return df, df_participants, subjects


def main():
    parser = get_parser()
    args = parser.parse_args()
    path_SC_1 = args.path_SC_1
    path_SC_2 = args.path_SC_2
    global name_SC_1
    name_SC_1 = args.name_SC_1
    global name_SC_2
    name_SC_2 = args.name_SC_2
    global ref
    ref = args.ref_subject
    vertlevels = args.vertlevel
    path_out_figures = os.path.join(args.path_out, 'figures')
    path_out_csv = os.path.join(args.path_out, 'csv')
    if args.angle:
        filename_angle = args.angle
    else:
        filename_angle = None
    # If the output folder directory is not     present, then create it.
    if not os.path.exists(path_out_figures):
        os.makedirs(path_out_figures)
    if not os.path.exists(path_out_csv):
        os.makedirs(path_out_csv)
    
    # Dump log file there
    log = os.path.join(args.path_out, FNAME_LOG)
    if os.path.exists(log):
        os.remove(log)
    fh = logging.FileHandler(os.path.join(log))
    logging.root.addHandler(fh)

    # Get exclude subjects list
    if args.exclude:
        exclude_list = args.exclude
    else: 
        exclude_list = None

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


    # Read csv files and create dataframe:
    df_1, df_participants, subjects = read_csv_files(path_SC_1, args.participant_file, exclude_list=exclude)
    df_1['method'] = name_SC_1
    df_2, df_participants, subjects = read_csv_files(path_SC_2, args.participant_file, exclude_list=exclude)
    df_2['method'] = name_SC_2

    # Find common participants
    common_participants = set(df_1['participant_id']).intersection(df_2['participant_id'])

    # Filter both DataFrames to include only common participants
    df1_filtered = df_1[df_1['participant_id'].isin(common_participants)]
    df2_filtered = df_2[df_2['participant_id'].isin(common_participants)]

    # Concatenate the filtered DataFrames
    df = pd.concat([df1_filtered, df2_filtered], ignore_index=True)
    
    path_out = path_out_figures
    # Compute descriptive statistics (mean and std age, weight, height)
    df = df.dropna(axis=1, how='all')
    subjects_after_dropping = df['Filename'].unique()

    # Print number of subjects
    logger.info(f'Number of subjects after dropping NaN: {str(len(subjects_after_dropping))}')
    logger.info(f'Dropped subjects: {str(list(set(list(subjects)) - set(list(subjects_after_dropping))))}\n')

    # Keep only VertLevel from C1 to Th1 (slices 705 to 969)
    df = df[df['Slice (I->S)'] <= 948]
    df = df[df['Slice (I->S)'] >= 722]
    # Multiply solidity by 100 to get percentage (sct_process_segmentation computes solidity in the interval 0-1)
    df['MEAN(solidity)'] = df['MEAN(solidity)'] * 100

    # Create plots
    create_lineplot(df, 'method', path_out)        # across all subjects
    create_lineplot(df.loc[df['method'] == name_SC_1], hue=None, path_out=path_out, filename='lineplot_' + name_SC_1 + '.png')        # across all subjects
    create_lineplot(df.loc[df['method'] == name_SC_2], hue=None, path_out=path_out, filename='lineplot_' + name_SC_2 + '.png')        # across all subjects
    # Create line plot of STD:
    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    axs = axes.ravel()
    for index, metric in enumerate(METRICS):
        hue = 'method'
        color_palette = PALETTE['method']
         # Compute std for each slice across all subjects
        sd_df = (df.groupby(["Slice (I->S)", hue])[metric]
            .std()
            .reset_index(name="sd"))
        sns.lineplot(
            ax=axs[index],
            data=sd_df,
            x="Slice (I->S)",
            y="sd",
            hue=hue,
            palette=color_palette,
            linewidth=2,
            legend=False
        )
        axs[index].invert_xaxis()
        axs[index].set_ylabel(f'STD of {METRIC_TO_TITLE[metric]}', fontsize=LABELS_FONT_SIZE)
        axs[index].set_xlabel(f'Axial Slice #', fontsize=LABELS_FONT_SIZE)
        axs[index].tick_params(axis='both', which='major', labelsize=TICKS_FONT_SIZE)
        axs[index].spines['right'].set_visible(False)
        axs[index].spines['left'].set_visible(False)
        axs[index].spines['top'].set_visible(False)
        axs[index].spines['bottom'].set_visible(True)

        # Get indices of slices corresponding vertebral levels
        vert, ind_vert, ind_vert_mid = get_vert_indices(df)
        # Insert a vertical line for each intervertebral disc
        for idx, x in enumerate(ind_vert[1:-1]):
            axs[index].axvline(df.loc[x, 'Slice (I->S)'], color='black', linestyle='--', alpha=0.5, zorder=0)
        ymin, ymax = axs[index].get_ylim()
        # Insert a text label for each vertebral level
        for idx, x in enumerate(ind_vert_mid, 0):
            # Deal with T1 label (C8 -> T1)
            if vert[x] > 8:
                level = 'T' + str(vert[x] - 8)
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)

            else:
                level = 'C' + str(vert[x])
                axs[index].text(df.loc[ind_vert_mid[idx], 'Slice (I->S)'], ymin, level, horizontalalignment='center',
                                verticalalignment='bottom', color='black', fontsize=TICKS_FONT_SIZE)

            # Add only horizontal grid lines
        axs[index].yaxis.grid(True)
        # Move grid to background (i.e. behind other elements)
        axs[index].set_axisbelow(True)

    fig.savefig(os.path.join(path_out, 'lineplot_std.png'), dpi=300, bbox_inches='tight')
    # Compute statistics between methods
    logger.info(f'\nComparing methods...')
    compare_metrics_across_methods(df)
    df_scaling = compute_scaling_factor(df)
    create_lineplot(df_scaling, hue=None, path_out=path_out, set_axis=False, filename='lineplot_scaling_factor.png')

    if filename_angle:
        import statsmodels.formula.api as smf
        from scipy.stats import pearsonr, spearmanr
        import statsmodels.api as sm
        df_angle = pd.read_csv(filename_angle)
        logger.info(f'\nComparing methods with angle as covariate...')
        # print mean and std of neck angle        
        logger.info(f'Neck angle: mean = {df_angle["Angle (degrees)"].mean()}, std = {df_angle["Angle (degrees)"].std()}')
        # Merge angle to df1_filtered: assign each subject's angle to all their slices
        angle_dict = dict(zip(df_angle['Subject'], df_angle['Angle (degrees)']))
        df['angle_deg'] = df['participant_id'].map(angle_dict)
        df = df.dropna(subset=['participant_id','angle_deg','method','MEAN(area)','Slice (I->S)']).reset_index(drop=True)        # Merge angle to df1_filtered
        m = smf.mixedlm('Q("MEAN(area)") ~ angle_deg + Q("Slice (I->S)") + C(method)', df, groups=df["participant_id"])
        res = m.fit()

        df["method"] = pd.Categorical(df["method"], categories=["discs", "spinal_level"])
        m = smf.mixedlm(
            'Q("MEAN(area)") ~ angle_deg * C(method) + Q("Slice (I->S)")',
            df,
            groups=df["participant_id"]
        )
        res = m.fit()
        logger.info(res.summary())

        df_roi = df[df['Slice (I->S)'].between(725, 800)].copy()
        #df_roi = df[df['Slice (I->S)'].between(800, 900)].copy()
        df_mean = (df_roi
                        .groupby(["participant_id", "method"], as_index=False)
                        .agg(mean_area=("MEAN(area)", "mean"),
                                angle_deg=("angle_deg", "first")))
        # print mean and std of angle_deg for each method:
        for meth in sorted(df_mean["method"].unique()):
            d = df_mean[df_mean["method"] == meth].dropna(subset=["angle_deg"])
            mean_angle = d["mean_area"].mean()
            std_angle = d["mean_area"].std()
            logger.info(f"Method: {meth} - mean_area mean = {mean_angle:.2f}, std = {std_angle:.2f}")

        for meth in sorted(df_mean["method"].unique()):
            d = df_mean[df_mean["method"] == meth].dropna(subset=["mean_area","angle_deg"])
            # check for normality of mean_area and angle_deg using Shapiro-Wilk test
            stat, p_value_mean_area = stats.shapiro(d["mean_area"])
            stat, p_value_angle_deg = stats.shapiro(d["angle_deg"])
            logger.info(f"\nMethod: {meth} - Shapiro-Wilk test for mean_area: p-value{format_pvalue(p_value_mean_area)}, Shapiro-Wilk test for angle_deg: p-value{format_pvalue(p_value_angle_deg)}")
            r, p = pearsonr(d["angle_deg"], d["mean_area"])
            rho, p_s = spearmanr(d["angle_deg"], d["mean_area"])
            logger.info(f"\nMethod: {meth}  (n={len(d)})")
            logger.info(f"  Pearson r={r:.3f}, p={p:.4g}")
            logger.info(f"  Spearman rho={rho:.3f}, p={p_s:.4g}")

        plt.figure(figsize=(6, 4), dpi=150)
        ax = sns.scatterplot(
            data=df_mean,
            x="angle_deg",
            y="mean_area",
            hue="method",
            markers="o",
            palette=color_palette,
            style="method",
            s=70
        )

        # optional: add per-method linear fit lines
        sns.regplot(
            data=df_mean[df_mean["method"] == "discs"],
            x="angle_deg", y="mean_area",
            scatter=False, ax=ax, color='blue', line_kws={"lw": 2}
        )
        sns.regplot(
            data=df_mean[df_mean["method"] == "spinal_level"],
            x="angle_deg", y="mean_area",
            scatter=False, ax=ax, color='red', line_kws={"lw": 2}
        )

        ax.set_xlabel("Neck angle [degrees]")
        ax.set_ylabel("Mean CSA [mm²]")
        ax.set_title("Mean CSA vs neck angle (per method)")
        ax.legend(title="Method", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        # Add text of spearman correlation permethod with pvalue formatted with format_pvalue function
        for meth in sorted(df_mean["method"].unique()):
            d = df_mean[df_mean["method"] == meth].dropna(subset=["mean_area","angle_deg"])
            rho, p_s = spearmanr(d["angle_deg"], d["mean_area"])
            text = f"rho = {rho:.3f}, p{format_pvalue(p_s)}"
            if meth == "discs":
                plt.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=9, verticalalignment='top', color='blue')
            else:
                plt.text(0.05, 0.90, text, transform=ax.transAxes, fontsize=9, verticalalignment='top', color='red')
        # save if you want
        plt.savefig(os.path.join(path_out,"angle_vs_mean_area_by_method.png"), bbox_inches="tight", dpi=300)
        logger.info(f'Created: {os.path.join(path_out,"angle_vs_mean_area_by_method.png")}\n')
plt.show()    
if __name__ == '__main__':
    main()
