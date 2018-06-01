rxnFile = open('/tmp/rxns.tsv','r')
for line in rxnFile:
    lps = line.split('\t')
    print(lps[5].strip(), lps[17].strip(), lps[18].strip())

