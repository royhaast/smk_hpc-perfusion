# lookup tables for setting wb structure:
hemi_to_structure = {'L': 'CORTEX_LEFT', 'Lflip': 'CORTEX_LEFT', 'R': 'CORTEX_RIGHT'}
surf_to_secondary_type = {'midthickness': 'MIDTHICKNESS', 'inner': 'PIAL', 'outer': 'GRAY_WHITE'}

# Generate subject-specific warps to go from a template unfolded
# space to the subject's native space (and vice versa). Does this 
# for inner, midthickness and outer levels
def construct_warp_cmd (wildcards, input):
    singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

    set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
    set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

    cmd = f"{set_matlab_lic} {set_java_home} {singularity_cmd} "\
            f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
            f"create_warps_old('{input.io_folder}','{input.io_folder}'); "\
            f"extrapolate_warp_unfold2native('{input.io_folder}','{input.io_folder}');\"" 
    return cmd   

rule create_unfold_to_native_warps:
    input:
        dseg = rules.autotop.output,
        io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop',
    output:
        native2unfold = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_native2unfold.nii',
        unfold2native = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_unfold2native.nii'
    params:
        matlab_cmd = construct_warp_cmd
    group: 'surfaces'        
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
    group: 'surfaces'     
    threads: 8
    resources:
        mem_mb = 32000    
    run:
        for in_unfolded, out_native in zip(input.unfolded, output.native):
            shell("wb_command -surface-apply-warpfield {in_unfolded} {input.warp} {out_native}")
            shell("wb_command -set-structure {out_native} CORTEX_RIGHT -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS")

# This uses the original surfaces from Jordan's code and crops them
# to the same resolution as the surfaces generated using the newer
# (Ali's) code
def construct_crop_cmd (wildcards, input):
    singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

    set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
    set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

    cmd = f"{set_matlab_lic} {set_java_home} {singularity_cmd} "\
            f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
            f"crop_SurfMap('{input.io_folder}','{input.io_folder}');\"" 
    return cmd   

rule crop_original_surfaces:
    input:
        dseg = rules.autotop.output,
        io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop'
    output:
        expand('results/autotop-dev/work/autotop/sub-{{subject}}/sub-{{subject}}_hemi-{{H}}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-{depth}.surf.gii', depth=['1','2','3','4'])
    params:
        matlab_cmd = construct_crop_cmd
    group: 'surfaces' 
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "{params.matlab_cmd}"

# Generate midthickness surface using connectome workbench instead
# of using the laplace coordinates
rule generate_midthickness_surf:
    input:
        inner = rules.crop_original_surfaces.output[0], #'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-1.surf.gii',
        outer = rules.crop_original_surfaces.output[3] #'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-4.surf.gii'
    output: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midSurf_depth-mid.surf.gii',
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.H]
    group: 'surfaces'     
    singularity: config['singularity_connectomewb'] 
    shell:
        "wb_command -set-structure {input.inner} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type PIAL && "
        "wb_command -set-structure {input.outer} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE && "
        "wb_command -surface-cortex-layer {input.outer} {input.inner} 0.5 {output}  && "
        "wb_command -set-structure {output} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS"

# Compute average warp across all subjects, per hemisphere to 
# generate group average hippocampal surfaces
def construct_avg_warp_cmd (wildcards,input,output):
    cmd = "fslmaths {}".format(input[0])
    for in_volume in input[1:]:
        cmd = "{0} -add {1}".format(cmd, in_volume)
    cmd = "{0} -div {1} {2}".format(cmd, len(input), output)
    return cmd    

rule create_avg_warp:
    input: expand('results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{{H}}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/Warp_unfold2native.nii', subject=subjects)
    output: 'results/average_hippocampus/Warp_unfold2native_hemi-{H}.nii.gz'
    params:
        average_cmd = construct_avg_warp_cmd
    group: 'surfaces'     
    singularity: config['singularity_prepdwi']
    shell:
        "{params.average_cmd}"   

rule transform_unfolded_to_avg_native:
    input:
        warp = rules.create_avg_warp.output,
        unfolded = 'resources/{surface}.unfolded.surf.gii'
    output: 'results/average_hippocampus/{surface}_hemi-{H}.nativecrop.surf.gii'
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.H],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surface]
    group: 'surfaces'     
    singularity: config['singularity_connectomewb']        
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "wb_command -surface-apply-warpfield {input.unfolded} {input.warp} {output} && "
        "wb_command -set-structure {output} {params.structure_type} -surface-type ANATOMICAL -surface-secondary-type {params.secondary_type}"