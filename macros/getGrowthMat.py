#!../venvpy/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/29     **
#**  last modified: 2018/05/29     **
#************************************

import argparse
import subprocess
import collections
import pandas

def main():
    args = options()
    verbose = args.verbose
    wte = args.exclude.split(' ')
    df = pandas.read_csv(args.biolog, sep='\t')
    print(df.shape)
    dfe = df[df['Well'].isin(wte)]
    print('Removing the following wells from growth matrix:', dfe)
    df = df[~df['Well'].isin(wte)]
    print(df.shape)
    gw = list(df.loc[df[args.id]>0, 'Well'])
    ngw = list(df.loc[df[args.id]==0, 'Well'])
    odf = pandas.read_csv(args.compounds, sep='\t')
    odf['Growth'] = -1
    odf.loc[odf['Well'].isin(gw), 'Growth'] = 1
    odf.loc[odf['Well'].isin(ngw), 'Growth'] = 0
    odf.to_csv(args.outpath+'growthMatrix_'+args.id+'.csv', index=False)

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-c', '--compounds', help='name of source compound', default='../rhizobiumRoot491/media/singleNMedia/ncompounds.tsv')
    parser.add_argument('-b', '--biolog', help='name of biolog tsv file', default='../rhizobiumRoot491/data/biolog_summary.tsv')
    parser.add_argument('-o', '--outpath', help='path for output file', default='../rhizobiumRoot491/data/')
    parser.add_argument('-i', '--id', help='id for the organism', default='Root491')
    parser.add_argument('-e', '--exclude', help='wells to be excluded from the growth matrix (to be used for validation)', default='A2 A5 A12 B6 B10')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
