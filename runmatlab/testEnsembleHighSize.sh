mkdir -p ../outputs/root491/hs/
nohup matlab -nodesktop -nosplash -nodisplay -r "ensembleHighSize;exit" > /tmp/hs.log 2>/tmp/errhs.out &
