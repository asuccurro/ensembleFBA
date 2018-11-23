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

shortened = {'C4-C9-D2-D6-D7-D10-E6-E7-E10-G4-G6-G7-G10-G12': 'not_found_and_G12',
             'A2-A5-A12-B6-B10-C4-C9-D2-D6-D7-D10-E6-E7-E10-G4-G6-G7-G10-G12': '5N_and_nf_and_G12'}

def main():
    args = options()
    verbose = args.verbose
    excl = args.exclude
    selc = args.select
    classid = args.classid
    orgid = args.id
    outpath = args.outpath
    biolog = args.biolog
    compounds = args.compounds

    idf, df = getGrowthMatrix(biolog, compounds, outpath, orgid, classid, excl, selc, verbose)
    return idf, df

def getGrowthMatrix(biolog, compounds, outpath, orgid, classid, excl, selc, verbose=False):

    wte = excl.split(' ')
    wts = selc.split(' ')
    
    oname = ''
    if len(wte[0]) > 0:
        estr = '-'.join(wte)
        if len(estr) > 30:
            if estr in shortened:
                oname = '_exclude_%s' % shortened[estr]
            else:
                print('ACHTUNG! String %s is too long, MATLAB will complain' % estr)
                oname = '_exclude_XYZ'
        else:
            oname = '_exclude_%s' % estr
    if len(wts[0]) > 0:
        oname = '_select_%s' % ('-'.join(wts))
    idf = pandas.read_csv(biolog, sep='\t')
    if classid != 'NA':
        if verbose:
            print('Selecting all ', classid, ' for growth matrix')
        wts = list(idf[idf["Class"] == classid]["Well"])
        oname = oname+'_select_%s' % (classid.replace(' ', '-'))
    df  = idf[~idf['Well'].isin(wte)]
    if verbose:
        print(idf.shape)
        print('Removing the following wells from growth matrix:', idf[idf['Well'].isin(wte)])
        print(df.shape)
    if len(wts[0]) > 0:
        df = df[df['Well'].isin(wts)]
        if verbose:
            print('Selecting the following wells from growth matrix:', idf[idf['Well'].isin(wts)])
            print(df.shape)
    gw = list(df.loc[df[orgid]>0, 'Well'])
    ngw = list(df.loc[df[orgid]==0, 'Well'])
    odf = pandas.read_csv(compounds, sep='\t')
    odf['Growth'] = -1
    odf.loc[odf['Well'].isin(gw), 'Growth'] = 1
    odf.loc[odf['Well'].isin(ngw), 'Growth'] = 0
    odf.to_csv(outpath+'growthMatrix_'+orgid+oname+'.csv', index=False)
    return idf, df

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-c', '--compounds', help='name of source compound', default='../rhizobiumRoot491/media/singleNMedia/ncompounds.tsv')
    parser.add_argument('-b', '--biolog', help='name of biolog tsv file', default='../rhizobiumRoot491/data/biolog_summary.tsv')
    parser.add_argument('-o', '--outpath', help='path for output file', default='../rhizobiumRoot491/data/')
    parser.add_argument('-i', '--id', help='id for the organism', default='Root491')
    parser.add_argument('-e', '--exclude', help='wells to be excluded from the growth matrix (to be used for validation)', default='A2 A5 A12 B6 B10')
    parser.add_argument('-s', '--select', help='wells to be selected for the growth matrix', default='')
    parser.add_argument('-l', '--classid', help='class of wells to be selected for the growth matrix', default='NA')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
