rule warp_perf_to_corobl_crop:
    input:
        nii = lambda wildcards: 'results/perfusion/sub-{s}/sub-{s}_{f}'.format(s=wildcards.subject, f=config['perfusion_data'][wildcards.asl_parameter]),
        init = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = join(config['hippunfold_dir'],config['template'])
    output: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_{hemi}.nii.gz'
    group: 'map_perfusion'
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm} -t {input.init}"

rule lr_flip_perf:
    input: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_L.nii.gz'
    output: 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_Lflip.nii.gz'
    group: 'map_perfusion'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

rule sample_perf_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_{H}.nii.gz',
        midthickness = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{asl_parameter}_{H}.nativecrop.shape.gii'
    group: 'map_perfusion'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -trilinear"

rule sample_perf_hippocampus_myelin_style:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_{asl_parameter}_{H}.nii.gz',
        ribbon = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/ribbon.nii.gz',
        midthickness = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii',
        thickness = 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.nativecrop.shape.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{asl_parameter}_myelin-style_{H}.nativecrop.shape.gii'
    group: 'map_perfusion'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -myelin-style {input.ribbon} {input.thickness} 0.6"

     