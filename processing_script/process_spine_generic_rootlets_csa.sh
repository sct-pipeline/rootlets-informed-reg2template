#!/bin/bash
#
# Process data of different neck positions (extension, flexion and straight).
#
# Usage:
#   ./process_data.sh <SUBJECT>
# 
# Manual segmentations and labels (discs, PMJ, nerve rootlets) should be located under:
# PATH_DATA/derivatives/labels/SUBJECT/anat/
#
# Authors: Sandrine Bédard

set -x
# Immediately exit if error
set -e -o pipefail

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT

# Retrieve input params
SUBJECT=$1

# Save script path
PATH_DERIVATIVES="${PATH_DATA}/derivatives/labels"

# get starting time:
start=`date +%s`

# FUNCTIONS
# ==============================================================================

# NOTE: manual disc labels should go from C1-C2 to C7-T1.
label_if_does_not_exist(){
  local file="$1"
  local file_seg="$2"
  # Update global variable with segmentation file name
  # Remove space other in filename 
  #suffix="_space-other"
  FILELABEL="${file//$suffix/}_label-discs_dlabel"
  FILELABELMANUAL="${PATH_DERIVATIVES}/${SUBJECT}/anat/${FILELABEL}.nii.gz"
  echo "Looking for manual label: $FILELABELMANUAL"
  if [[ -e $FILELABELMANUAL ]]; then
    echo "Found! Using manual labels."
    rsync -avzh $FILELABELMANUAL ${FILELABEL}.nii.gz
    # Generate labeled segmentation from manual disc labels
    sct_image -i ${file}.nii.gz -set-sform-to-qform
    sct_image -i ${file_seg}.nii.gz -set-sform-to-qform
    sct_image -i ${FILELABEL}.nii.gz -set-sform-to-qform
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -discfile ${FILELABEL}.nii.gz -c t2 -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    echo "Not found. Proceeding with automatic labeling."
    # Generate labeled segmentation
    sct_label_vertebrae -i ${file}.nii.gz -s ${file_seg}.nii.gz -c t2 -qc ${PATH_QC} -qc-subject ${SUBJECT}
  fi
}

# Check if manual segmentation already exists. If it does, copy it locally. If it does not, perform segmentation.
segment_if_does_not_exist(){
  local file="$1"
  folder_contrast="anat"

  # Update global variable with segmentation file name
  FILESEG="${file}_seg"
  FILESEGMANUAL="${PATH_DATA}/derivatives/labels_softseg_bin/${SUBJECT}/${folder_contrast}/${file}_desc-softseg_label-SC_seg.nii.gz"
  echo
  echo "Looking for manual segmentation: $FILESEGMANUAL"
  if [[ -e $FILESEGMANUAL ]]; then
    echo "Found! Using manual segmentation."
    rsync -avzh $FILESEGMANUAL ${FILESEG}.nii.gz
    sct_qc -i ${file}.nii.gz -s ${FILESEG}.nii.gz -p sct_deepseg_sc -qc ${PATH_QC} -qc-subject ${SUBJECT}
  else
    echo "Not found. Proceeding with automatic segmentation."
    # Segment spinal cord
    sct_deepseg -i ${file}.nii.gz spinalcord -qc ${PATH_QC} -largest 1 -qc-subject ${SUBJECT} -o ${FILESEG}.nii.gz
  fi
}

# Check if manual rootlets already exists. If it does, copy it locally. If it does not, perform segmentation.
segment_rootlets_if_does_not_exist(){
  local file="$1"
  # Update global variable with segmentation file name
  FILELABEL="${file}_rootlets"
  set +e
  FILELABELMANUAL="$(ls ${PATH_DERIVATIVES}/${SUBJECT}/anat/*rootlets*.nii.gz)"
  FILELABELMANUAL="$FILELABEL"  # TODO: earase this line if want manual rootelts to be used
  set -e -o pipefail
  echo "Looking for manual label: $FILELABELMANUAL"
  if [[ -e $FILELABELMANUAL ]]; then
    echo "Found! Using manual rootlets labels."
    rsync -avzh $FILELABELMANUAL ${file}_rootlets.nii.gz
  else
    echo "Not found. Proceeding with automatic labeling."
    # Segment rootlets
    sct_deepseg -i ${file}.nii.gz rootlets -qc ${PATH_QC} -o ${FILELABEL}.nii.gz

  fi
}

# SCRIPT STARTS HERE
# ==============================================================================
# Display useful info for the log, such as SCT version, RAM and CPU cores available
sct_check_dependencies -short

# Go to folder where data will be copied and processed
cd ${PATH_DATA_PROCESSED}
# Copy source images
rsync -Ravzh ${PATH_DATA}/./${SUBJECT}/anat/${SUBJECT}_*T2w.* .
# Go to anat folder where all structural data are located
cd ${SUBJECT}/anat/
file_t2="${SUBJECT}_T2w"

# 1. Segment spinal cord
segment_if_does_not_exist $file_t2

# 3. Segment rootlets
segment_rootlets_if_does_not_exist $file_t2
file_t2_rootlets=$FILELABEL

# 4. Compute CSA based on rootlets
mkdir -p ${PATH_RESULTS}/rootlets_PAM50
sct_process_segmentation -i ${file_t2}_seg.nii.gz -rootlet ${file_t2_rootlets}.nii.gz -perslice 1 -normalize-PAM50 1 -o ${PATH_RESULTS}/rootlets_PAM50/${file_t2}_PAM50.csv -qc ${PATH_QC} -qc-subject ${SUBJECT}

mkdir -p ${PATH_RESULTS}/rootlets_perlevel
sct_process_segmentation -i ${file_t2}_seg.nii.gz -rootlet ${file_t2_rootlets}.nii.gz -perslice 1 -normalize-PAM50 1 -o ${PATH_RESULTS}/rootlets_perlevel/${file_t2}_PAM50.csv -qc ${PATH_QC} -qc-subject ${SUBJECT}


# Display useful info for the log
# ===============================
end=`date +%s`
runtime_final=$((end-start))
echo
echo "~~~"
echo "SCT version: `sct_version`"
echo "Ran on:      `uname -nsr`"
echo "Duration:    $(($runtime_final / 3600))hrs $((($runtime_final / 60) % 60))min $(($runtime_final % 60))sec"
echo "~~~"
