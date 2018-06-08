for c in exc5N inc5N; do
	mkdir -p ../outputs/root491/${c}
	for s in 0 1; do
		for n in 21 81; do
			for t in 0; do
				fnm=ensemble_${c}_${s}_${n}_${t}
				cp ../macros/00_ensemble.m ${fnm}.m
				sed -i "s/XXXTEST/${t}/g" ${fnm}.m
				sed -i "s/XXXSTOC/${s}/g" ${fnm}.m
				sed -i "s/XXXSIZE/${n}/g" ${fnm}.m
				sed -i "s/XXXCOND/${c}/g" ${fnm}.m
				nohup matlab -nodesktop -nosplash -nodisplay -r "${fnm};exit" > /tmp/${fnm}.log 2>/tmp/err_${fnm}.out &
			done
		done
	done
done
