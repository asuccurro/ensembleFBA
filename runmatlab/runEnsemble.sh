matlab -nodesktop -nosplash -nodisplay -r "ensemble;exit" > /tmp/ensemble.log 2>/tmp/err.out 
mkdir -p ../outputs/ensemble/
mv Ensemble* ../outputs/ensemble/
