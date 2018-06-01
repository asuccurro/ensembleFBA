for c in cpd17041 cpd17042 cpd17043 cpd12370; do
echo $c
grep $c ../rhizobiumRoot491/db/reactions.tsv > /tmp/rxns.tsv
python checkRxns.py
done

