# Calculate different morphological features using the surface
# geometry. Includes gyrification (i.e., vertex area), 
# curvature and thickness
rule calculate_gyrification_from_surface:
    input: 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_gyrification_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb'] 
    group: 'morphology' 
    shell:
        "wb_command -surface-vertex-areas {input} {output}"    

rule calculate_curvature_from_surface:
    input: 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_curvature_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-curvature {input} -mean {output}"     

rule calculate_thickness_from_surface:
    input:
        inner = 'results/surface_warps/sub-{subject}/{H}/inner.native.surf.gii',
        outer = 'results/surface_warps/sub-{subject}/{H}/outer.native.surf.gii',
    output: 'results/surface_maps/sub-{subject}/sub-{subject}_thickness_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'       
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"  

# Calculate distance from subject-specific surfaces to group
# template. Does this for midthickness only.
rule calculate_distance_to_avg:
    input:
        avg = 'results/average_hippocampus/midthickness_hemi-{H}.native.surf.gii',
        subject = 'results/surface_warps/sub-{subject}/{H}/midthickness.native.surf.gii'
    output:
        distance = 'results/surface_maps/sub-{subject}/sub-{subject}_surfdist_{H}.native.shape.gii',
        displacement = 'results/surface_maps/sub-{subject}/sub-{subject}_surfdisp_{H}.native.shape.gii'
    singularity: config['singularity_connectomewb']
    group: 'morphology'              
    shell:    
        "wb_command -surface-to-surface-3d-distance {input.avg} {input.subject} {output.distance} -vectors {output.displacement} && "
        "wb_command -metric-math 'metric' {output.distance} -fixnan 0 -var metric {output.distance} && "
        "wb_command -metric-math 'metric' {output.displacement} -fixnan 0 -var metric {output.displacement}"