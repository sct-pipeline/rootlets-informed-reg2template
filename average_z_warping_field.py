#!/usr/bin/env python
# -*- coding: utf-8
# TODO
#
# For usage, type: python average_z_warping_field.py.py -h

# Authors: Sandrine Bédard

import argparse
import numpy as np
import nibabel as nib
import os
import numpy as np


def get_parser():
    parser = argparse.ArgumentParser(
        description="Averages Z warping field for each slice, does not count zeros in the mean.")
    parser.add_argument('-i', required=True, type=str,
                        help="Warping field in Z direction only. Split before.")
    parser.add_argument('-o', required=False, type=str,
                        default='warp_Z_mean.nii.gz',
                        help="Output warping field filename.")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    src = nib.load(args.i)
    src_np = np.array(src.get_fdata())
    # Change zeros to nan to ignore during mean
    src_np[src_np == 0] = np.nan
    src_np_mean = np.zeros_like(src_np)
   # print(src_np_mean)
    for slices in range(src_np.shape[2]):
        if np.isnan(src_np[:,:,slices]).all():
            mean_np = [[0]]
        else:
            mean_np = np.nanmean(src_np[:,:,slices],keepdims=True)
        src_np_mean[:,:,slices] = np.broadcast_to(mean_np, src_np[:,:,slices].shape)
    nii_mean_warp= nib.Nifti1Image(src_np_mean, src.affine)
    fname_out_levels = args.o
    nib.save(nii_mean_warp, fname_out_levels)


if __name__ == '__main__':
    main()
