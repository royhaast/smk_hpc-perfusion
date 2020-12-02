# Prepare TOF data for mapping on unfolded hippocampus
rule warp_tof_to_corobl_crop:
    input:
        nii = 'results/tof/sub-{subject}/sub-{subject}_acq-TOF_run-avg_angio_reg2TSE_Warped.nii.gz',
        init = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = join(config['hippunfold_dir'],config['template'])
    output: 'results/maps/sub-{subject}/sub-{subject}_TOF_{hemi}.nii.gz'
    group: 'map_tof'
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm} -t {input.init}"

rule lr_flip_tof:
    input: 'results/maps/sub-{subject}/sub-{subject}_TOF_L.nii.gz'
    output: 'results/maps/sub-{subject}/sub-{subject}_TOF_Lflip.nii.gz'
    group: 'map_tof'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

rule sample_tof_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_TOF_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        inner = rules.generate_midthickness_surf.input.inner,
        midthickness = rules.generate_midthickness_surf.output,
        outer = rules.generate_midthickness_surf.input.outer
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_TOF_{H}.nativecrop.shape.gii'
    group: 'map_tof'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -ribbon-constrained {input.outer} {input.inner} -volume-roi {input.ribbon}"

# Prepare MP2RAGE data for reconstruction of vasculature
rule tof_n4_correction:
    input: 'results/maps/sub-{subject}/sub-{subject}_TOF_{H}.nii.gz'
    output: 'results/vasculature/sub-{subject}/{H}/sub-{subject}_TOF_{H}_n4.nii.gz'
    group: 'map_vasculature'
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "N4BiasFieldCorrection -d 3 -i {input} -o {output}"    

rule tof_stedi_filter:
    input: rules.tof_n4_correction.output
    output: 'results/vasculature/sub-{subject}/{H}/sub-{subject}_TOF_{H}_n4_stedi.nii.gz'
    group: 'map_vasculature'
    #singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "bash scripts/tof_filter.sh {input}"      

rule tof_extract_vessels:
    input: rules.tof_stedi_filter.output
    output: 'results/vasculature/sub-{subject}/{H}/sub-{subject}_TOF_{H}_n4_stedi_centerdia_clean.nii.gz'
    group: 'map_vasculature'
    singularity: config['singularity_braincharter']
    threads: 8
    resources:
        mem_mb = 32000        
    shell:
        "bash scripts/tof_extract_vessels.sh $(realpath {input}) 'nii.gz' TOF"

""" 
Manual work needed from here for 3D reconstruction of the vasculature using the 
output files of the 'tof_extract_vessels' or 'tof_n4_correction' rule and MeVisLab. 
A template MeVisLab file (*.mlab) can be found in the 'resources' directory. After, 
generating the vessel masks (i.e., '') run the 'tof_extract_vessels' rule

"""
# Extract features from reconstructed vasculature
# rule tof_binarize_segmentation:
#     input: 'results/vasculature/sub-{subject}/sub-{subject}_TOF_n4_vessel_seg_1p0.nii.gz'
#     output: 'results/vasculature/sub-{subject}/sub-{subject}_TOF_n4_vessel_seg_1p0_bin.nii.gz'
#     group: 'map_tof'
#     singularity: config['singularity_prepdwi']
#     shell:
#         "fslmaths {input} -bin {output}"       

rule tof_extract_diameters:
    input: 'results/vasculature/sub-{subject}/{H}/sub-{subject}_TOF_{H}_n4_vessel_seg_0p7.nii.gz'
    output: 'results/vasculature/sub-{subject}/{H}/sub-{subject}_TOF_{H}_n4_vessel_seg_1p0_diameters_clean_dil.nii.gz'
    group: 'map_vasculature'
    singularity: config['singularity_braincharter']   
    threads: 8
    resources:
        mem_mb = 32000     
    shell:
        "bash scripts/tof_extract_diameters.sh $(realpath {input}) 'nii.gz'"