#! /bin/bash

## Master script for MPRAGE correction
# 1. Mirror BIDS folder structure for MP2RAGE and Sa2RAGE images
# 2. Run main functionality

s=$1

# Mirror BIDS folder
gradcorrect=$2 #"/scratch/rhaast/HPC_perfusion/data/bids"
mp2rage_correction="/scratch/rhaast/HPC_perfusion/results/mp2rage_correction"

# if [ ! -d $mp2rage_correction ] ; then
#     mkdir -p $mp2rage_correction
# fi

pushd $gradcorrect
#cp `ls -p | grep -v /` $mp2rage_correction/sub-${s}

mkdir -p $mp2rage_correction/sub-${s}/anat && cp -nr anat/*MP2RAGE* $mp2rage_correction/sub-${s}/anat/
mkdir -p $mp2rage_correction/sub-${s}/fmap && cp -nr fmap/*SA2RAGE* $mp2rage_correction/sub-${s}/fmap/

popd

# Run MP2RAGE correction
script_home=`dirname $0`

bash $script_home/hpc_MP2RAGE_preprocessing.sh $s

gzip $mp2rage_correction/sub-${s}/anat/*.nii
gzip $mp2rage_correction/sub-${s}/fmap/*.nii

# # Run gradient distortion correction
# singularity_gradcorrect="/project/6007967/akhanf/singularity/bids-apps/khanlab_gradcorrect_v0.0.3a.sif"
# gradcorrect="/scratch/rhaast/HPC_perfusion/data/gradcorrect"
# coeff="/scratch/rhaast/MSTRCHT_BIDSify/scripts/coeff_SC72CD.grad"

# singularity run $SINGULARITY_OPTS $singularity_gradcorrect $mp2rage_correction $gradcorrect participant --grad_coeff_file $coeff --participant_label $s
