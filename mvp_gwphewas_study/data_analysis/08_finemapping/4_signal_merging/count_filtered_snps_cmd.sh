#!/usr/bin/bash
path=$1
command_path=$2
out_path=$3

#Reset command path
>$command_path

for i in $path/*.filtered_out.txt; do
    echo "wc___-l___"$i"___>>___""$out_path" >> $command_path
done

