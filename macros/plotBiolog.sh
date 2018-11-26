source ../venvpy/bin/activate

f=ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1
for o in Root9 Root491 Root66D1; 
	do
	for c in _exclude_not_found_and_G12 _exclude_5N_and_nf_and_G12; 
	do
	p="../outputs/$o$c/"
	echo $p
	python makeBiologFigure.py -f $f --iopath $p -M
	python makeBiologFigure.py -U -f $f --iopath $p -M
	echo ""
	done
done

deactivate
