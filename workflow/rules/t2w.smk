rule combine_lr_manual_seg:
    input:
        lh = join(config['data']['manual_segs'],'sub-{subject}_hemi-L_hipp.nii.gz'),
        rh = join(config['data']['manual_segs'],'sub-{subject}_hemi-R_hipp.nii.gz'),
    output: join(config['data']['manual_segs'],'sub-{subject}_hemi-LR_hipp.nii.gz'),
    group: 't2w'
    run: 
        import numpy as np
        import nibabel as nib

        L  = nib.load(input.lh)
        R  = nib.load(input.rh)
        LR = L.get_fdata()+R.get_fdata()

        img = nib.Nifti1Image(LR.astype(np.int8), header=L.header, affine=L.affine)
        nib.save(img, output[0])

rule relabel_lr_manual_seg:
    input: rules.combine_lr_manual_seg.output
    output: join(config['data']['manual_segs'],'sub-{subject}_hemi-LR_hipp_relabeled.nii.gz')
    group: 't2w'    
    run:
        import numpy as np
        import nibabel as nib

        LR = nib.load(input[0])
        LR_labels = LR.get_fdata()
        LR_relabeled = np.zeros(LR_labels.shape)

        relabeling = [[1, 1],[2, 2],[20, 3],
                      [21, 4],[22, 5],[23, 6]]

        for labels in relabeling:
            LR_relabeled[LR_labels==labels[0]] = labels[1]

        img = nib.Nifti1Image(LR_relabeled.astype(np.int8), header=LR.header, affine=LR.affine)
        nib.save(img, output[0])

rule get_xfm_manual_seg_to_t2w:
    input:
        flo = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0.nii.gz'),
        ref = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-TSE_run-01_T2w.nii.gz')
    output:
        xfm_ras = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-ras_xfm.txt'),
        xfm_itk = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_from-dseg_to-refT2w_type-itk_xfm.txt'),
        warped = join(config['data']['manual_segs'],'sub-{subject}_acq-TSE_0p3_template0_space-refT2w.nii.gz')
    group: 't2w'        
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000     
    shell:
        'reg_aladin -flo {input.flo} -ref {input.ref} -res {output.warped} -aff {output.xfm_ras} -rigOnly -nac && '
        'c3d_affine_tool {output.xfm_ras} -oitk {output.xfm_itk}'

rule apply_xfm_manual_seg_to_t2w:
    input:
        nii = rules.relabel_lr_manual_seg.output,
        ref = rules.get_xfm_manual_seg_to_t2w.input.flo,
        xfm = rules.get_xfm_manual_seg_to_t2w.output.xfm_itk
    output: join(config['data']['manual_segs'],'sub-{subject}_desc-manualseg.nii.gz')
    group: 't2w'    
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output} -r {input.ref} -t {input.xfm}'     

rule smooth_labels:
    input: rules.apply_xfm_manual_seg_to_t2w.output
    output: join(config['data']['manual_segs'],'sub-{subject}_desc-manualseg_dseg.nii.gz')
    group: 't2w'    
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "c3d {input} -split -foreach -smooth 0.3mm -endfor -merge -o {output}"

rule autotop:
    input:
        bids = '/scratch/rhaast/MSTRCHT_BIDSify/gradcorrect',
        dseg = rules.smooth_labels.output,
        t1w = rules.copy_depadded.output.nii
    output: 'results/autotop-dev/results/sub-{subject}/anat/sub-{subject}_space-T1w_desc-subfields_modality-segT2w_template-CITI168_dseg.nii.gz'
    params:
        out_dir = 'results/autotop-dev',
        autotop = config['autotop'],
        autotop_opts = '--derivatives /scratch/rhaast/HPC_perfusion/data/manual_segs --modality segT2w --filter_T1w acquisition=MP2RAGEdepadded',
        snakemake_opts = '--cores 8 --use-singularity -p &> logs/autotop/sub-{subject}.log'
    group: 'autotop'        
    threads: 8
    resources:
        time = 120,
        mem_mb = 32000
    #log: ''
    shell:
        "{params.autotop} {input.bids} {params.out_dir} participant --participant_label {wildcards.subject} {params.autotop_opts} {params.snakemake_opts}"

def get_matlab_cmd (wildcards, input):
    singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

    set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
    set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

    cmd = f"{set_matlab_lic} {set_java_home} {singularity_cmd} "\
            f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
            f"create_warps('{input.io_folder}','{input.io_folder}')\""
    return cmd   

rule create_unfold_to_native_warps:
    input:
        dseg = rules.autotop.output,
        io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop',
    output:
        native2unfold = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_native2unfold.nii',
        unfold2native = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_unfold2native.nii'
    params:
        matlab_cmd = get_matlab_cmd
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "{params.matlab_cmd}"

rule transform_unfolded_to_native:
    input:
        warp = rules.create_unfold_to_native_warps.output.unfold2native,
        unfolded = expand('resources/{surface}.unfolded.surf.gii',surface=['inner','midthickness','outer'])
    output:
        native = expand('results/autotop-dev/work/autotop/sub-{{subject}}/sub-{{subject}}_hemi-{{H}}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/{surface}.nativecrop.surf.gii',surface=['inner','midthickness','outer'])
    #singularity: config['singularity_connectomewb']        
    threads: 8
    resources:
        mem_mb = 32000    
    run:
        for in_unfolded, out_native in zip(input.unfolded, output.native):
            shell("wb_command -surface-apply-warpfield {in_unfolded} {input.warp} {out_native}")
            shell("wb_command -set-structure {out_native} CORTEX_RIGHT -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS")

rule warp_t2w_to_corobl_crop:
    input:
        nii = join(config['data']['gradcorrect'],'anat/sub-{subject}_acq-TSE_rec-avg.nii.gz'),
        init = rules.get_xfm_manual_seg_to_t2w.output.xfm_itk,
        xfm = 'results/autotop-dev/work/sub-{subject}/anat/sub-{subject}_desc-affine_from-T2w_to-CITI168corobl_type-itk_xfm.txt',
        ref = join(config['hippunfold_dir'],config['template'])
    output: 'results/maps/sub-{subject}/sub-{subject}_T2w_{hemi}.nii.gz'
    group: 'map_t2w'
    singularity: config['singularity_prepdwi']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.nii} -o {output} -r {input.ref}  -t {input.xfm} -t {input.init}"

rule lr_flip_t2w:
    input: 'results/maps/sub-{subject}/sub-{subject}_T2w_L.nii.gz'
    output: 'results/maps/sub-{subject}/sub-{subject}_T2w_Lflip.nii.gz'
    group: 'map_t2w'
    singularity: config['singularity_prepdwi']
    shell:
        "c3d {input} -flip x -o  {output}"

rule sample_t2w_hippocampus:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_T2w_{H}.nii.gz',
        midthickness = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_T2w_{H}.nativecrop.shape.gii'
    group: 'map_t2w'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -trilinear"

rule extract_gm_ribbon:
    input: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/coords-IO.nii.gz'
    output: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/ribbon.nii.gz'
    singularity: config['singularity_prepdwi']
    shell:
        "fslmaths {input} -bin {output}"

rule sample_t2w_myelin_style:
    input:
        nii = 'results/maps/sub-{subject}/sub-{subject}_T2w_{H}.nii.gz',
        ribbon = rules.extract_gm_ribbon.output,
        midthickness = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii',
        thickness = 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.nativecrop.shape.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_T2w_myelin-style_{H}.nativecrop.shape.gii'
    group: 'map_t2w'
    singularity: config['singularity_connectomewb']
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        "wb_command -volume-to-surface-mapping {input.nii} {input.midthickness} {output} -myelin-style {input.ribbon} {input.thickness} 0.6"