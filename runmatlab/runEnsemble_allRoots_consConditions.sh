ngc=10
gc=21
for o in Root9 Root491 Root66D1; do
    for c in "" _exclude_A2-A5-A12-B6-B10; do
	mkdir -p ../outputs/${o}${c}
	for s in 1; do
	    for n in 50; do
		for t in 0; do
		    fnmx=ensemble_${o}${c}_${n}_${gc}_${ngc}_${s}_${t}
		    ## Because matlab does not like - passed through command line (interpreted as operation!)
		    fnm=${fnmx//-/_}
		    cp ../macros/template_ensemble.m ${fnm}.m
		    sed -i "s/XXXTEST/${t}/g" ${fnm}.m
		    sed -i "s/XXXSTOC/${s}/g" ${fnm}.m
		    sed -i "s/XXXSIZE/${n}/g" ${fnm}.m
		    sed -i "s/XXXCOND/${c}/g" ${fnm}.m
		    sed -i "s/XXXORG/${o}/g" ${fnm}.m
		    sed -i "s/XXXNGC/${ngc}/g" ${fnm}.m
		    sed -i "s/XXXGC/${gc}/g" ${fnm}.m
		    nohup matlab -nodesktop -nosplash -nodisplay -r "${fnm};exit" > /tmp/${fnm}.log 2>/tmp/err_${fnm}.out </dev/null &
		done
	    done
	done
    done
done
