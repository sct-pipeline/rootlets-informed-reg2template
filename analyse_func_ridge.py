#!/usr/bin/env python
# -*- coding: utf-8

# fslstats thresh_zstat1_med.nii.gz -V | cut -d " " -f1 >> voxels_med.txt
# fslstats thresh_zstat1_high.nii.gz -M >> zscore_high.txt
# fslstats thresh_zstat1_med.nii.gz -H 15 3.1 9 >> hitogram_med_15_3.1-9.txt

import os
import argparse
import nibabel as nib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i1",
                        required=True,
                        type=str,
                        help="First zstats file")
    parser.add_argument("-i2",
                        required=True,
                        type=str,
                        help="Second zstats file")
    parser.add_argument('-o', required=False, type=str,
                        help="Output filename")
    return parser


def load_voxel_values_by_z(file_path):
    """Load voxel values and group them by Z-axis slices."""
    img = nib.load(file_path)
    data = img.get_fdata()

    if len(data.shape) != 3:
        raise ValueError("Expected a 3D fMRI activation map (x, y, z). Got shape: {}".format(data.shape))

    z_slices = data.shape[2]  # Number of Z slices
    voxel_dict = {}

    for z in range(z_slices):
        slice_data = data[:, :, z].flatten()
        slice_data = slice_data[~np.isnan(slice_data)]  # Remove NaNs
        slice_data = slice_data[slice_data != 0]  # Remove zero values
        voxel_dict[z] = slice_data

    return voxel_dict


def create_ridge_plot(nii_file1, nii_file2, output_file):
    """Generate and save a ridge plot with Z-axis as the Y-axis for two input files."""
    sns.set(style="white", context="talk")

    voxel_dict1 = load_voxel_values_by_z(nii_file1)
    voxel_dict2 = load_voxel_values_by_z(nii_file2)

    data1 = [(z, v, "Rootlets-reg") for z, voxels in voxel_dict1.items() for v in voxels]
    data2 = [(z, v, "Discs-reg") for z, voxels in voxel_dict2.items() for v in voxels]

    df = pd.DataFrame(data1 + data2, columns=["Z-Slice", "Voxel Intensity", "File"])

    # Plot ridge density along the Z-axis for both files
    plt.figure(figsize=(3, 12))  # Adjust figure size for better fit
    sns.kdeplot(data=df[df["File"] == "Rootlets-reg"], y="Z-Slice", weights="Voxel Intensity", fill=True, alpha=0.7, bw_adjust=0.5, color="red", label="Rootlets-reg")
    sns.kdeplot(data=df[df["File"] == "Discs-reg"], y="Z-Slice", weights="Voxel Intensity", fill=True, alpha=0.7, bw_adjust=0.5, color="royalblue", label="Discs-reg")
    plt.legend(loc="upper center")
    # Formatting
    plt.ylabel("")
    plt.xlabel("Active Voxel Density")
    plt.xlim(0, 0.035)  # Fix x-axis limits to 0.035
    plt.gca().xaxis.set_tick_params(which='both', direction='in', length=6, width=1, color='black')  # Add tick lines with customization
    plt.ylim(58, 234)
    #plt.yticks(np.arange(58, 245, 10))  # Adjust y-axis ticks for better readability
    plt.gca().yaxis.set_ticks([])  # Remove y-axis ticks
    plt.gca().xaxis.set_tick_params(which='both', direction='in', length=6, width=1, color='black')  # Add tick lines for x-axis
    sns.despine(left=True, bottom=False)  # Remove the box around the plot
    plt.xticks(fontsize=16)  # Adjust x-axis label font size
    plt.xticks(np.arange(0, 0.036, 0.01))  # Adjust x-axis ticks for better readability
    #plt.yticks(fontsize=16)  # Adjust y-axis label font size
    plt.tight_layout()  # Ensure everything fits within the figure

    # Save the figure
    plt.savefig(output_file, dpi=300)
    print(f"Ridge plot saved as: {output_file}")


def main():
    args = get_parser().parse_args()
    print(f"Input file 1: {args.i1}")
    print(f"Input file 2: {args.i2}")
    create_ridge_plot(args.i1, args.i2, os.path.join(args.o, 'ridge_plot.png'))


if __name__ == "__main__":
    main()
