#!usr/bin/bash
script=$1
script_dir=`readlink -f $script`
script_dir=$(dirname "${script_dir}")
script=`basename $script .m`

in=$2

out=$3
out_dir=`dirname $out`

in_cat12=$out_dir/mp2rage.nii

if (file $in | grep -q compressed ) ; then
    gunzip -k -c $in > $in_cat12
else
    cp $in $in_cat12    
fi

module load matlab ants

pushd $script_dir
    matlab -nosplash -nodisplay -nodesktop -r "${script}('${in_cat12}','${out_dir}'); quit"
popd

cat_out=$out_dir/mri/p0mp2rage
ThresholdImage 3 $cat_out.nii ${out} 1.2 4 1 0
