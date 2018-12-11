source ../venvpy/bin/activate

for c in M; do
    # U give empty!
    for t in C; do
	python getEssentialGenesStats.py -D -${t} -${c}
    done
done

deactivate
