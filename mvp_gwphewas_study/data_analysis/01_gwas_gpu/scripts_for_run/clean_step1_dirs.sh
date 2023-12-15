#!/usr/bin/bash
path=".//output/HARE_ANC_Run"
group=$1
for i in $path/*/$group/step1; do
    if [ ! -f $i/*.out.gpu.rda ]; then
	#echo "NO RDA"
	rm $i/phewas.ld.maf.*.gpu*
        echo $i;
    elif [ ! -s $i/*.out.gpu.rda ]; then
	#echo "RDA: 0byte"
        rm $i/phewas.ld.maf.*.gpu*
        echo $i;
    elif [ ! -f $i/*.out.gpu.varianceRatio.txt ]; then
	#echo "NO VAR"
        rm $i/phewas.ld.maf.*.gpu*
	echo $i;
    elif [ ! -s $i/*.out.gpu.varianceRatio.txt ]; then
	#echo "VAR: 0byte"
        rm $i/phewas.ld.maf.*.gpu*
        echo $i;
    elif [ ! -f $i/*gpu_30markers.SAIGE.results.txt ]; then 
	#echo "NO MARKER"
        rm $i/phewas.ld.maf.*.gpu*
	echo $i;
    elif [ `grep closed $i/*.gpu.log | wc -l` == 0 ]; then
	#echo "NO closed"
        rm $i/phewas.ld.maf.*.gpu*
	echo $i;
    fi
done

