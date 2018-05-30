import pandas

df = pandas.read_csv('../data/ModelSEEDdata/compounds.tsv', sep='\t')
df = df[['id','name']]
df.to_csv('simple-cpd-DB.tsv', sep='\t', index=False)

df = pandas.read_csv('../data/ModelSEEDdata/reactions.tsv', sep='\t')
df = df[['id','name']]
df.to_csv('simple-rxn-DB.tsv', sep='\t', index=False)
