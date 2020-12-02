# Calculate different morphological features using the surface
# geometry. Includes gyrification (i.e., vertex area), 
# curvature and thickness
rule calculate_gyrification_from_surface:
    input: rules.generate_midthickness_surf.output
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_gyrification_{H}.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb'] 
    group: 'morphology' 
    shell:
        "wb_command -surface-vertex-areas {input} {output}"    

rule calculate_curvature_from_surface:
    input: rules.generate_midthickness_surf.output
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_curvature_{H}.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-curvature {input} -mean {output}"     

rule calculate_thickness_from_surface:
    input:
        inner = rules.generate_midthickness_surf.input.inner,
        outer = rules.generate_midthickness_surf.input.outer
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"  

# Calculate distance from subject-specific surfaces to group
# template. Does this for midthickness only.
rule calculate_distance_to_avg:
    input:
        avg = 'results/average_hippocampus/midthickness_hemi-{H}.nativecrop.surf.gii',
        subject = 'results/autotop-dev/work/autotop/sub-{subject}/sub-{subject}_hemi-{H}_space-CITI168corobl_desc-cropped_modality-segT2w_autotop/midthickness.nativecrop.surf.gii'
    output:
        distance = 'results/surface_maps/sub-{subject}/sub-{subject}_surfdist_{H}.nativecrop.shape.gii',
        displacement = 'results/surface_maps/sub-{subject}/sub-{subject}_surfdisp_{H}.nativecrop.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'              
    shell:    
        "wb_command -surface-to-surface-3d-distance {input.avg} {input.subject} {output.distance} -vectors {output.displacement} && "
        "wb_command -metric-math 'metric' {output.distance} -fixnan 0 -var metric {output.distance} && "
        "wb_command -metric-math 'metric' {output.displacement} -fixnan 0 -var metric {output.displacement}"