def get_morphology_cmd (wildcards, input, output):
    singularity_cmd = f"singularity exec -B {config['autotop_matlab']}:/src -e {config['singularity_autotop']}" 

    set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
    set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"

    out_thickness = os.path.realpath(output.thickness)
    out_gyrification = os.path.realpath(output.gyrification)

    cmd = f"{set_matlab_lic} {set_java_home} {singularity_cmd} "\
            f"{config['matlab_bin']} -batch \"addpath(genpath('{config['autotop_matlab']}')); "\
            f"map_gm_morphology_gifti('{input.io_folder}','{out_thickness}','{out_gyrification}')\""
    return cmd   

rule get_morphology_maps:
    input:
        dseg = rules.autotop.output,
        io_folder = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop',
    output:
        thickness = 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.nativecrop.shape.gii',
        gyrification = 'results/surface_maps/sub-{subject}/sub-{subject}_gyrification_{H}.nativecrop.shape.gii'
    params:
        matlab_cmd = get_morphology_cmd
    threads: 8
    resources:
        mem_mb = 32000    
    shell:
        "SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "{params.matlab_cmd}"

rule calculate_curvature:
    input: 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_curvature_{H}.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']  
    shell:
        "wb_command -surface-curvature {input} -mean {output}"     

# def construct_merge_cmd (wildcards,input,output):
#     cmd = "wb_command -metric-merge {}".format(output)
#     for in_metric in input:
#         cmd = "{0} -metric {1} -column 1".format(cmd, in_metric)
#     return cmd

# rule merge_t2w_maps:
#     input: expand('results/surface_maps/sub-{subject}/sub-{subject}_T2w_{{H}}.nativecrop.shape.gii', subject=subjects)
#     output: 'results/surface_maps/group/sub-group_T2w_{H}.nativecrop.shape.gii'
#     params:
#         merge_cmd = construct_merge_cmd
#     group: 'map_t2w'
#     singularity: config['singularity_connectomewb']  
#     shell:
#         "{params.merge_cmd}"    

# rule average_t2w_maps:
#     input: rules.merge_t2w_maps.output 
#     output:
#         avg = 'results/surface_maps/group/sub-group_T2w_{H}_avg.nativecrop.shape.gii',
#         std = 'results/surface_maps/group/sub-group_T2w_{H}_std.nativecrop.shape.gii',
#         cov = 'results/surface_maps/group/sub-group_T2w_{H}_cov.nativecrop.shape.gii'
#     group: 'map_t2w'
#     singularity: config['singularity_connectomewb']  
#     shell:
#         "wb_command -metric-reduce {input} MEAN {output.avg} -only-numeric &&"
#         "wb_command -metric-reduce {input} STDEV {output.std} -only-numeric &&"
#         "wb_command -metric-reduce {input} COV {output.cov} -only-numeric"