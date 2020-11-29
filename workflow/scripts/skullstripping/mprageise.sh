#! /bin/bash

####
# Written by Sriranga Kashyap <kashyap.sriranga@gmail.com>
# Maastricht, 2020-06
####

export AFNI_NIFTI_TYPE_WARN=NO
export AFNI_ENVIRON_WARNINGS=NO

# Simple formatting
bold=$(tput bold)
normal=$(tput sgr0)

# Help
function Help() {
    cat <<HELP

Usage:

$(basename $0) ${bold}-i${normal} INV2 image ${bold}-u${normal} UNI image ${bold}-o${normal} out folder ${bold}-r${normal} re-bias

--------------------------------------------------------------------------------
Input arguments:

    -i: MP2RAGE INV2 magnitude image ( e.g. /path/to/source/inv2_mag.nii.gz )

    -u: MP2RAGE UNI image ( e.g. /path/to/source/uni.nii.gz )

    -o: Output folder ( e.g. /path/to/output )

    -r: Reintroduce bias-field ( default = 0 )

--------------------------------------------------------------------------------

Example:

$(basename $0) -i /data/inv2.nii.gz -u /data/uni.nii.gz

--------------------------------------------------------------------------------
Script was created by: Sriranga Kashyap (06-2020), kashyap.sriranga@gmail.com
--------------------------------------------------------------------------------

HELP
    exit 1
}

# Check for flag
if [[ "$1" == "-h" || $# -eq 0 ]]; then
    Help >&2
fi

# Get some info
Aversion=$(afni -ver)
runDate=$(echo $(date))


# Parse input arguments
while getopts "h:i:r:o:u:" OPT; do
    case $OPT in
        h) #help
            Help
            exit 0
        ;;
        i) # INV2 image
            inv2_image=$OPTARG
        ;;
        r) # Rebias boolean
            re_bias=$OPTARG
        ;;
        o) # Output folder
            out_folder=$OPTARG
        ;;        
        u) # UNI image
            uni_image=$OPTARG
        ;;
        \?) # report error
            echo "$HELP" >&2
            exit 1
        ;;
    esac
done

if [ ! -z "$re_bias" ]; then
    
    echo "++ UNI will ${bold} NOT ${normal} be bias-corrected."

    re_bias=1
    
else
    
    echo "++ UNI will be bias-corrected."
    
    re_bias=0
    
fi

# Check extension
if [ ${inv2_image##*.} = gz ]; then
    inv2_basename=$(basename $inv2_image .nii.gz)
    uni_basename=$(basename $uni_image .nii.gz)
else
    inv2_basename=$(basename $inv2_image .nii)
    uni_basename=$(basename $uni_image .nii)
fi

# Check if bias field should be removed or not
if [ "$re_bias" = 0 ]; then
    
    echo "++ Removing bias-field."
    
    # Unifize INV2
    3dUnifize \
    -quiet \
    -overwrite \
    -prefix $out_folder/${inv2_basename}_bfc.nii.gz \
    $inv2_image
    
else
    
    echo "++ Reintroducing bias-field."
    
    3dcopy $inv2_image $out_folder/${inv2_basename}_bfc.nii.gz
    
fi

# Intensity Normalise INV2 between 0-1
int_max=$(3dinfo -dmaxus $inv2_image)
int_min=$(3dinfo -dminus $inv2_image)

3dcalc \
-overwrite \
-a $out_folder/${inv2_basename}_bfc.nii.gz \
-expr "( a - $int_min ) / ( $int_max - $int_min )" \
-prefix $out_folder/${inv2_basename}_intnorm.nii.gz

echo "++ MPRAGEising the UNI image."

# Multiply intensity normalised INV2 to UNI image
3dcalc \
-overwrite \
-a $uni_image \
-b $out_folder/${inv2_basename}_intnorm.nii.gz \
-expr "a * b" \
-prefix $out_folder/${uni_basename}_unbiased_clean.nii.gz

# Remove intermediate files
rm $out_folder/${inv2_basename}_bfc.nii.gz $out_folder/${inv2_basename}_intnorm.nii.gz

echo "++ Done."

