nohup matlab -nodesktop -nosplash -nodisplay -r "gapfill;exit" > /tmp/gapfill.log 2>/tmp/err.out &
mkdir -p ../outputs/gapfill_sequential
mv GapFill_Sequence* ../outputs/gapfill_sequential/
