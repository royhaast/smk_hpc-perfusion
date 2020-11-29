def construct_merge_cmd (wildcards,input,output):
    cmd = "wb_command -metric-merge {}".format(output)
    for in_metric in input:
        cmd = "{0} -metric {1} -column 1".format(cmd, in_metric)
    return cmd

rule merge_group_metric_gii:
    input: expand('results/surface_maps/sub-{subject}/sub-{subject}_{{parameter}}_{{H}}.nativecrop.shape.gii', subject=subjects)
    output: 'results/surface_maps/group/sub-group_{parameter}_{H}.nativecrop.shape.gii'
    params:
        merge_cmd = construct_merge_cmd
    singularity: config['singularity_connectomewb']  
    shell:
        "{params.merge_cmd}"    

rule reduce_group_metric_gii:
    input: rules.merge_group_metric_gii.output 
    output:
        avg = 'results/surface_maps/group/sub-group_{parameter}_{H}_avg.nativecrop.shape.gii',
        std = 'results/surface_maps/group/sub-group_{parameter}_{H}_std.nativecrop.shape.gii',
        cov = 'results/surface_maps/group/sub-group_{parameter}_{H}_cov.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']  
    shell:
        "wb_command -metric-reduce {input} MEAN {output.avg} -only-numeric &&"
        "wb_command -metric-reduce {input} STDEV {output.std} -only-numeric &&"
        "wb_command -metric-reduce {input} COV {output.cov} -only-numeric "   

# concatenate subject's L and R data
rule merge_subject_lr_metric_gii:
    input: expand('results/surface_maps/sub-{{subject}}/sub-{{subject}}_{{parameter}}_{H}.nativecrop.shape.gii', H=['Lflip','R'])
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{parameter}_LR.nativecrop.shape.gii'
    params:
        merge_cmd = construct_merge_cmd
    singularity: config['singularity_connectomewb']  
    shell:
        "{params.merge_cmd}"

# average across subject's L and R data
rule reduce_subject_lr_metric_gii:
    input: rules.merge_subject_lr_metric_gii.output 
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_{parameter}_LR_avg.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']  
    shell:
        "wb_command -metric-reduce {input} MEAN {output} -only-numeric"

# concatenate all subject's L and R data
rule merge_group_lr_metric_gii:
    input: expand('results/surface_maps/sub-{subject}/sub-{subject}_{{parameter}}_{H}.nativecrop.shape.gii', subject=subjects, H=['Lflip','R'])
    output: 'results/surface_maps/group/sub-group_{parameter}_LR.nativecrop.shape.gii'
    params:
        merge_cmd = construct_merge_cmd
    singularity: config['singularity_connectomewb']  
    shell:
        "{params.merge_cmd}"      

# average across all subject's L and R data
rule reduce_group_lr_metric_gii:
    input: rules.merge_group_lr_metric_gii.output 
    output:
        avg = 'results/surface_maps/group/sub-group_{parameter}_LR_avg.nativecrop.shape.gii',
        std = 'results/surface_maps/group/sub-group_{parameter}_LR_std.nativecrop.shape.gii',
        cov = 'results/surface_maps/group/sub-group_{parameter}_LR_cov.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']  
    shell:
        "wb_command -metric-reduce {input} MEAN {output.avg} -only-numeric &&"
        "wb_command -metric-reduce {input} STDEV {output.std} -only-numeric &&"
        "wb_command -metric-reduce {input} COV {output.cov} -only-numeric"   