# Prepare ASL data for mapping on unfolded hippocampus
rule warp_perf_to_corobl_crop:
    input:
        nii = lambda wildcards: 'results/perfusion_sdc/sub-{s}/sub-{s}_{f}.nii.gz'.format(s=wildcards.subject, f=config['perfusion_data'][wildcards.asl_parameter]),
        xfm2anat = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_M0-to-MP2RAGE_BBR_fsl2greedy2ants_Apply1st.txt',
        xfm2tse = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_MP2RAGE-to-TSE_final_greedy2ants_Apply2nd.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = join(config['hippunfold_dir'],config['template'])
    output: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_{hemi}.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']   
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse} -t {input.xfm2anat}"

rule lr_flip_perf:
    input: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_L.nii.gz'
    output: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_Lflip.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

rule sample_perf_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{asl_parameter}_{H}.native.shape.gii'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_connectomewb']   
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"

# Prepare individual ASL runs for mapping on unfolded hippocampus
rule warp_run_mean_to_corobl_crop:
    input:
        nii = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_acq-ASL_run-{run}_moCorr_sDC_PWI_Tmean.nii.gz',
        xfm2anat = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_M0-to-MP2RAGE_BBR_fsl2greedy2ants_Apply1st.txt',
        xfm2tse = 'results/perfusion_sdc/sub-{subject}/sub-{subject}_MP2RAGE-to-TSE_final_greedy2ants_Apply2nd.txt',
        xfm2ref = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm2crop = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = join(config['hippunfold_dir'],config['template'])
    output: 'results/perfusion_preprocessing/sub-{subject}/sub-{subject}_run-{run}_mean_{hemi}.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']   
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm2crop} -t {input.xfm2ref} -t {input.xfm2tse} -t {input.xfm2anat}"

rule lr_flip_run_mean:
    input: 'results/perfusion_preprocessing/sub-{subject}/sub-{subject}_run-{run}_mean_L.nii.gz'
    output: 'results/perfusion_preprocessing/sub-{subject}/sub-{subject}_run-{run}_mean_Lflip.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

rule merge_run_means:
    input: expand('results/perfusion_preprocessing/sub-{{subject}}/sub-{{subject}}_run-{run}_mean_{{H}}.nii.gz', run=[str(r).zfill(2) for r in range(1,9)])
    output: 'results/maps/sub-{subject}/sub-{subject}_PWI_means_{H}.nii.gz'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_prepdwi']
    shell:
        "fslmerge -t {output} {input}"    

rule sample_run_means_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_PWI_means_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        midthickness = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_PWI_means_{H}.native.shape.gii'
    group: 'map_perfusion_hpc'
    singularity: config['singularity_connectomewb']   
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"

# Transform ASL data to lores MP2RAGE for cortical analyses
rule warp_perf_to_anatomy:
    input:
        unpack(collect_input),
        nii = lambda wildcards: 'results/perfusion_clean/sub-{s}/sub-{s}_{f}.nii.gz'.format(s=wildcards.subject, f=config['perfusion_data'][wildcards.asl_parameter]),
        xfm = 'results/perfusion_clean/sub-{subject}/sub-{subject}_acq-MZeroScan_distCorr_template_reg2Anat_greedy2ants_apply1st.txt'
    output: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_space-MP2RAGE.nii.gz'
    group: 'map_perfusion_ctx'
    singularity: config['singularity_prepdwi']   
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.t1w} -t {input.xfm}"

rule sample_perf_cortex:
    input:
        nii = rules.warp_perf_to_anatomy.output,
        inner = 'results/hcp_mmp/sub-{subject}/{hemi}.white.32k_fs_LR.surf.gii',
        midthickness = 'results/hcp_mmp/sub-{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii',
        outer = 'results/hcp_mmp/sub-{subject}/{hemi}.pial.32k_fs_LR.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{asl_parameter}_{hemi}.32k_fs_LR.shape.gii'
    group: 'map_perfusion_ctx'
    singularity: config['singularity_connectomewb']
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.inner} {input.outer} -voxel-subdiv 5"
