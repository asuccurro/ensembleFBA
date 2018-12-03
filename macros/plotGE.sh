source ../venvpy/bin/activate

for c in A M; do
    # U give empty!
    for t in O C; do
	python getEssentialGenesStats.py -D -${t} -${c}
    done
done

deactivate
