#!/usr/bin/sh

# creates 'data.in' input file for a given radial distance
template () {

    file=$1
    rd=$2

    > $file

    echo "He2" >> $file
    echo "2.0, 2.0, 0.0" >> $file
    echo "${rd}, 0.5" >> $file
    echo "0, 4" >> $file
    echo "15,1.5 , 15,1.5 , 15,1.5 , 15,1.5 , 15,1.5 , 15,1.5 , 15,1.5 , 15,1.5 , 15,1.5" >> $file
    echo "100.0, 20.0, 8, 32, 6, 20" >> $file
    echo "1.0e-20, 1.0e-20, 1.0e-20" >> $file
    echo "0, 0" >> $file
    echo "1, 0" >> $file
    echo "0, 0" >> $file
    echo "2, 100, 1.0e-10" >> $file
}

# ensures latest version of 'main' executable is used
cp "../code/main" "./hf_main"

# absolute parameters
input_dir="$(pwd)"
output_dir="$(pwd)"

data="${input_dir}/data.in"
results="${output_dir}/hf_results.dat"
curve="${output_dir}/curve.dat"

# command-line parameters
rd_min="$1"
rd_max="$2"
n="$3"

# radial-distance step
step="$(awk "BEGIN{ printf \"%.10f\n\", ($rd_max - $rd_min)/$n }")"

# reset 'curve.dat'
> $curve

# loop over radial-distance values
ii=0
while [ $ii -lt $n ]
do
    printf "${ii} of ${n}"

    rd="$(awk "BEGIN{ printf \"%.10f\n\", $rd_min + ($ii * $step) }")"

    template ${data} ${rd}

    ./hf_main "${input_dir}" "${output_dir}" > /dev/null

    echo ${rd} $(sed '6q;d' ${results}) >> "${curve}"

    true $(( ii++ ))
done
