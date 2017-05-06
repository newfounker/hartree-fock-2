#!/usr/bin/sh

# creates 'input.dat' input file for a given radial distance
template () {

    file=$1
    rd=$2
    ltop=$3

    > $file

    echo "He2" >> $file
    echo "2.0, 2.0, 0.0" >> $file
    echo "${rd}, 0.5" >> $file
    echo "0, ${ltop}" >> $file
    # echo "30,1.5 , 30,1.5 , 25,1.5 , 25,1.5 , 20,1.5 , 20,1.5 , 15,1.5 , 15,1.5 , 10,1.5" >> $file
    echo "15,1.5 ,  15,1.5 ,  15,1.5 ,  15,1.5 ,  15,1.5 ,  15,1.5 ,  15,1.5 ,  15,1.5" >> $file
    echo "100.0, 20.0, 8, 32, 6, 20" >> $file
    echo "1.0e-20, 1.0e-20, 1.0e-20" >> $file
    echo "0, 0" >> $file
    echo "1, 0" >> $file
    echo "0, 0" >> $file
    echo "2, 100, 1.0e-10" >> $file
}

# absolute parameters
input="$(pwd)/input.dat"
output_hf="$(pwd)/hf_results.dat"
output_curve="$(pwd)/curve.dat"

# command-line parameters
l_max="$1"
rd_min="$2"
rd_max="$3"
n="$4"

# radial-distance step
step="$(awk "BEGIN{ printf \"%.10f\n\", ($rd_max - $rd_min)/$n }")"

# loop over radial-distance values
ii=0
while [ $ii -lt ${n} ]
do
    rd="$(awk "BEGIN{ printf \"%.10f\n\", $rd_min + ($ii * $step) }")"

    printf "\r(${ii} / ${n}) ~ ${rd}"

    template ${input} ${rd} ${l_max}

    ../code/main ${input} ${output_hf} ${output_curve} > /dev/null

    true $(( ii++ ))
done

printf "\n"

# # loop over radial-distance values
# ii=0
# while [ $ii -lt ${n} ]
# do
#     rd="$(awk "BEGIN{ printf \"%.10f\n\", $rd_min + ($ii * $step) }")"

#     printf "${rd} " >> ${curve}

#     # loop over angular-momentum values
#     jj=0
#     while ! [ $jj -gt ${l_max} ]
#     do
#         printf "\r(${ii} / ${n}) , (${jj} / ${l_max})"

#         template ${data} ${rd} ${jj}

#         ./hf_main "${input_dir}" "${output_dir}" > /dev/null

#         # assuming energy is 7-th line in hf_results.dat
#         printf "$(sed '7q;d' ${results})" >> ${curve}

#         true $(( jj++ ))
#         true $(( jj++ ))
#     done

#     printf "\n" >> ${curve}
#     printf "\n"

#     true $(( ii++ ))
# done
