#!/usr/bin/bash
path=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_dir
output_path=/gpfs/alpine/med112/proj-shared/results/PostGWAS_Analysis/finemaping/GIA/matrix_counts_dir_EUR
for i in $path/*.map; do
    name=`basename $i`
    if [ ! -f $output_path/$name.count ]; then
        echo "wc___-l___"$i"___>___""$output_path/$name.count"
    fi
done
for i in $path/*.ld; do
    name=`basename $i`
    if [ ! -f $output_path/$name.count ]; then
	echo "wc___-l___"$i"___>___""$output_path/$name.count"
    fi
done
