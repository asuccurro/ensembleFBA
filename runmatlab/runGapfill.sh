nohup matlab -nodesktop -nosplash -nodisplay -r "gapfill;exit" > /tmp/gapfill.log 2>/tmp/err.out &
mkdir -p ../outputs/gapfill_sequential
mv Gapfill* ../outputs/gapfill_sequential/
