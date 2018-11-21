ngc=10
gc=21
for o in Root9 Root491 Root66D1; do
    #for c in "" _exclude_A2-A5-A12-B6-B10; do
    for c in _exclude_not_found_and_G12 _exclude_A2-A5-A12-B6-B10; do
	mydir="${o}${c}"
	fnm=ensemble_50_size_${gc}_gcs_${ngc}_ngcs_stochasticWeights_1
	sfnmx=ensemble_${o}${c}_50_${gc}_${ngc}_1_0
	sfnm=${sfnmx//-/_}
	cp ../macros/00_analyseEnsemble.m ae_${sfnm}.m
	sed -i "s/XXXFNAME/${fnm}/g" ae_${sfnm}.m
	sed -i "s@XXXDNAME@${mydir}@g" ae_${sfnm}.m
	echo $mydir
	echo $sfnm
	echo $fnm
	matlab -nodesktop -nosplash -nodisplay -r "ae_${sfnm};exit" 
    done
done

	
