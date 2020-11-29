rule mp2rage_correction:
    input: config['data']['gradcorrect']
    output: 'results/mp2rage_correction/sub-{subject}/anat/sub-{subject}_acq-MP2RAGE_run-01_corrT1_clean.nii.gz'
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/mp2rage_correction/sub-{subject}.log'
    shell:
        "bash scripts/mp2rage_correction/wrapper.sh {wildcards.subject} {input}"     

# rule tar_mp2rage_correction:
#     input: 'results/mp2rage_correction/sub-{subject}'
#     output: 'results/mp2rage_correction/sub-{subject}_mp2rage_corrected.tar.gz'
#     shell:
#         "cd results/mp2rage_correction && tar -czf `basename {output}` sub-{wildcards.subject}"

rule mprageise:
    input: unpack(collect_input)
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_run-01_corrUNI_clean_unbiased_clean.nii.gz'
    group: 'preprocessing'
    singularity: config['singularity_fmriprep']    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/mprageise/sub-{subject}.log'
    shell:
        "bash scripts/skullstripping/mprageise.sh -i {input.inv2} -u {input.t1w} -o `dirname {output}`" 

rule get_pads_mask:
    input: rules.mprageise.output #unpack(collect_input)
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_run-01_corrUNI_clean_unbiased_clean_dePadded.nii.gz'
    group: 'preprocessing'
    singularity: config['singularity_fmriprep']    
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/padsoff/sub-{subject}.log'    
    shell:
        "bash scripts/skullstripping/PadsOff -i {input}"

rule copy_depadded:
    input:
        nii = rules.get_pads_mask.output,
        json = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGE_run-01_UNI.json')
    output:
        nii = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGEdepadded_run-01_T1w.nii.gz'),
        json = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-MP2RAGEdepadded_run-01_T1w.json')
    shell:
        "cp {input.nii} {output.nii} && cp {input.json} {output.json}"


rule skull_stripping:
    input: rules.get_pads_mask.output #unpack(collect_input)
    output: 'results/skullstripping/sub-{subject}/sub-{subject}_acq-MP2RAGE_mask.nii.gz'
    params:
        script = 'scripts/skullstripping/run_cat12.m'
    group: 'preprocessing'
    threads: 8
    resources:
        mem_mb = 32000
    log: 'logs/skullstripping/sub-{subject}.log'
    shell:
        "bash scripts/skullstripping/skullstrip.sh {params.script} {input} `realpath {output}` " #&> {log} {params.out_dir}  

rule apply_brain_mask:
    input: 
        t1w = rules.mprageise.output,
        brain_mask = rules.skull_stripping.output
    output: 'results/freesurfer/sub-{subject}/mri/orig/001.mgz'
    group: 'preprocessing'
    singularity: config['singularity_freesurfer']
    shell:
        "mri_mask {input.t1w} {input.brain_mask} {output}"    

rule freesurfer:
    input: rules.apply_brain_mask.output
    output: 'results/freesurfer/sub-{subject}/scripts/recon-all.done'
    params:
        sd = 'results/freesurfer'
    group: 'freesurfer'
    singularity: config['singularity_freesurfer']
    threads: 8
    resources:
        time = 300,
        mem_mb = 32000
    shell:
        "export SUBJECTS_DIR={params.sd} && "
        "recon-all -all -s sub-{wildcards.subject} -no-wsgcaatlas -notal-check -threads 8"

rule extract_freesurfer_hippocampi:
    input:
        freesurfer = rules.freesurfer.output,
        t2w = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-TSE_rec-avg.nii.gz')
    output:
        lh = 'results/freesurfer/sub-{subject}/lh.hippoAmygLabels-T1-T2w.v21.mgz',
        rh = 'results/freesurfer/sub-{subject}/rh.hippoAmygLabels-T1-T2w.v21.mgz'
    params:
        sd = '/scratch/rhaast/HPC_perfusion/results/freesurfer'    
    group: 'postfreesurfer'
    singularity: config['singularity_freesurfer']
    threads: 8
    resources:
        time = 300,
        mem_mb = 32000
    shell:
        "segmentHA_T2.sh sub-{wildcards.subject} {input.t2w} T2w 1 {params.sd}"   

rule warp_t1_to_corobl_crop:
    input:
        nii = lambda wildcards: 'results/t1/sub-{s}/sub-{s}_{f}'.format(s=wildcards.subject, f=config['mp2rage_data'][wildcards.mp2rage_parameter]),
        init = 'data/manual_segs/sub-{subject}/anat/sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt',
        xfm = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = join(config['hippunfold_dir'],config['template'])
    output: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{hemi}.nii.gz'
    group: 'map_t1' 
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm} -t {input.init}"  

rule lr_flip_t1:
    input: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_L.nii.gz',
    output: 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_Lflip.nii.gz',
    group: 'map_t1'     
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

rule sample_t1_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{H}.nii.gz',
        midthickness = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{H}.nativecrop.shape.gii'
    group: 'map_t1'     
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -trilinear"

rule sample_t1_myelin_style:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_{H}.nii.gz',
        ribbon = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/ribbon.nii.gz',
        midthickness = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii',
        thickness = 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.nativecrop.shape.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{mp2rage_parameter}_myelin-style_{H}.nativecrop.shape.gii'
    group: 'map_t1'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -myelin-style {input.ribbon} {input.thickness} 0.6"

rule calculate_myelin_map:
    input:
        t1w = 'results/surface_maps/sub-{subject}/sub-{subject}_T1w_{H}.nativecrop.shape.gii',
        t2w = 'results/surface_maps/sub-{subject}/sub-{subject}_T2w_{H}.nativecrop.shape.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_myelin_{H}.nativecrop.shape.gii'    
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -metric-math 'clamp((t1w/t2w),0,5)' {output} -var t1w {input.t1w} -var t2w {input.t2w} -fixnan 0"