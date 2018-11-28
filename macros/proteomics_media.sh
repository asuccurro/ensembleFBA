mkdir -p ../data/MPIRoots/proteomicsMedia/
odir=../data/MPIRoots/proteomicsMedia/
for c in cpd00013 cpd00073 cpd00023 cpd00039 cpd00054;
do 
    cp ../data/MPIRoots/minimalMediaNoN.tsv $odir/mm_${c}.tsv
    printf "$c\t-100\t5\t0.001\n" >> $odir/mm_${c}.tsv
done
