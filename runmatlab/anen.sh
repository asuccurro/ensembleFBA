c1="exc5N inc5N"
array1=( $c1 )
c2="26 29"
array2=( $c2 )
for ((i=0;i<2;i++)); do
    mydir="root491/${array1[$i]}"
    g="${array2[$i]}"
    #for s in 0 1; do
    for s in 1; do
	#for n in 21 81; do
	for n in 21; do
	    fnm=ensemble_${n}_size_${g}_gcs_11_ngcs_stochasticWeights_${s}
	    cp ../macros/00_analyseEnsemble.m ae_${fnm}.m
	    sed -i "s/XXXFNAME/${fnm}/g" ae_${fnm}.m
	    sed -i "s@XXXDNAME@${mydir}@g" ae_${fnm}.m
	    echo $mydir
	    echo $fnm
	    matlab -nodesktop -nosplash -nodisplay -r "ae_${fnm};exit" 
	done
    done
done
