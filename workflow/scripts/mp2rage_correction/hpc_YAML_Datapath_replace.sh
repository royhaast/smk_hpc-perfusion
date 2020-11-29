#! /bin/bash
for s in {01..11}; do

    ## Low-Res MP2RAGE
    # MP2RAGE path
    sed -i "s#/scratch/rhaast/HPC_perfusion/data/mp2rage_correction/sub-${s}/anat#/scratch/rhaast/HPC_perfusion/results/mp2rage_correction/sub-${s}/anat#g" sub-${s}_lowres_mp2rage_config.yml
    # Sa2RAGE
    sed -i "s#/scratch/rhaast/HPC_perfusion/data/mp2rage_correction/sub-${s}/fmap#/scratch/rhaast/HPC_perfusion/results/mp2rage_correction/sub-${s}/fmap#g" sub-${s}_lowres_mp2rage_config.yml

    for r in {1..2}; do
        ## Hi-Res MP2RAGE
        # MP2RAGE path
        sed -i "s#/scratch/rhaast/HPC_perfusion/data/mp2rage_correction/sub-${s}/anat#/scratch/rhaast/HPC_perfusion/results/mp2rage_correction/sub-${s}/anat#g" sub-${s}_hires_mp2rage_run${r}_config.yml
        # Sa2RAGE
        sed -i "s#/scratch/rhaast/HPC_perfusion/data/mp2rage_correction/sub-${s}/fmap#/scratch/rhaast/HPC_perfusion/results/mp2rage_correction/sub-${s}/fmap#g" sub-${s}_hires_mp2rage_run${r}_config.yml

    done

done
