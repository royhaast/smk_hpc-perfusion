#! /bin/bash

s=$1

module load matlab ants
mat_bin="matlab -nodisplay -nodesktop"

utils_path="'/scratch/rhaast/HPC_perfusion/scripts/mp2rage_utils'"
spm_path="'/scratch/rhaast/HPC_perfusion/scripts/skullstripping/spm12'"

lowres_config="'/scratch/rhaast/HPC_perfusion/scripts/mp2rage_correction/sub-${s}_lowres_mp2rage_config.yml'"
$mat_bin -r "addpath(genpath($utils_path)); addpath($spm_path); mp2rage_main($lowres_config); exit;"

hires1_config="'/scratch/rhaast/HPC_perfusion/scripts/mp2rage_correction/sub-${s}_hires_mp2rage_run1_config.yml'"
uni_file="/scratch/rhaast/HPC_perfusion/results/mp2rage_correction/sub-${s}/anat/sub-${s}_acq-hiresMP2RAGE_run-01_UNI.nii.gz"

echo "ImageMath 3 ${uni_file} RescaleImage ${uni_file} 0 4095"
ImageMath 3 $uni_file RescaleImage $uni_file 0 4095

echo $mat_bin -r "addpath(genpath($utils_path)); addpath($spm_path); mp2rage_main($hires1_config); exit;"
$mat_bin -r "addpath(genpath($utils_path)); addpath($spm_path); mp2rage_main($hires1_config); exit;"

hires2_config="'/scratch/rhaast/HPC_perfusion/scripts/mp2rage_correction/sub-${s}_hires_mp2rage_run2_config.yml'"
uni_file="/scratch/rhaast/HPC_perfusion/results/mp2rage_correction/sub-${s}/anat/sub-${s}_acq-hiresMP2RAGE_run-02_UNI.nii.gz"

echo "ImageMath 3 ${uni_file} RescaleImage ${uni_file} 0 4095"
ImageMath 3 $uni_file RescaleImage $uni_file 0 4095

echo $mat_bin -r "addpath(genpath($utils_path)); addpath($spm_path); mp2rage_main($hires2_config); exit;"
$mat_bin -r "addpath(genpath($utils_path)); addpath($spm_path); mp2rage_main($hires2_config); exit;"