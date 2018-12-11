orgs=( Root9 Root491 Root66D1 )
gbid=( ASE33_ ASD46_ ASE09_ )
cpds=( cpd00013 cpd00023 cpd00039 cpd00054 cpd00073 )
name=( ammo glut lysi seri urea)
for ((i=0;i<3;i++)); do
    for ((j=0;j<5;j++)); do
	o=${orgs[$i]}
	g=${gbid[$i]}
	c=${cpds[$j]}
	n=${name[$j]}
	fnm=ge_${o}_${c}
	cp ../macros/template_geneEssAn.m ${fnm}.m
	sed -i "s/XXXORG/${o}/g" ${fnm}.m
	sed -i "s/XXXCPD/${c}/g" ${fnm}.m
	sed -i "s/XXXGID/${g}/g" ${fnm}.m
	sed -i "s/XXXNAM/${n}/g" ${fnm}.m
	nohup matlab -nodesktop -nosplash -nodisplay -r "${fnm};exit" > /tmp/${fnm}.log 2>/tmp/err_${fnm}.out </dev/null &
    done
done
