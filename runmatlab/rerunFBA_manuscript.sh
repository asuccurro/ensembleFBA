ngc=11
gc=26
for o in Root9 Root491 Root66D1; do
    for c in _exclude_not_found_and_G12 _exclude_5N_and_nf_and_G12; do
	mydir="${o}${c}"
	fnm=ensemble_50_size_${gc}_gcs_${ngc}_ngcs_stochasticWeights_1
	sfnmx=ensemble_${o}${c}_50_${gc}_${ngc}_1_0
	sfnm=${sfnmx//-/_}
	echo $mydir
	echo $sfnm
	echo $fnm
	cp ../macros/template_rerunFBA.m bm_${sfnm}.m
	sed -i "s/XXXFNAME/${fnm}/g" bm_${sfnm}.m
	sed -i "s@XXXDNAME@${mydir}@g" bm_${sfnm}.m
	sed -i "s/XXXCPDLIST/'cpd00013'; 'cpd00073'; 'cpd00023'; 'cpd00039'; 'cpd00054'/g" bm_${sfnm}.m
	sed -i "s/XXXCPDNAME/proteomics/g" bm_${sfnm}.m
	matlab -nodesktop -nosplash -nodisplay -r "bm_${sfnm};exit" 
    done
done

