for g in $1/*.log
do
       echo $g	
	for c in `cat $g` 
	do
		grep $c ../rhizobiumRoot491/data/growthMatrix_Root491.csv
	done
done
