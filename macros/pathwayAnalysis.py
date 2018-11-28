#!/usr/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/11/20     **
#**  last modified: 2018/11/20     **
#************************************

import argparse
import pandas
import numpy as np
import string
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns

colors={'L-Lysine': '#8DA0CB',
        'L-Glutamic-Acid': '#FC8D62',
        'L-Serine': '#E78AC3',
        'Urea': '#E78AC3'}

def mergeColumns(dfdic, coln, cols, thr=0):
    ndf = pandas.DataFrame(dfdic[cols[0]][coln])
    ndf.columns = [cols[0]]
    for c in cols[1:]:
        ndf[c] = dfdic[c][coln]
    #print(ndf.shape)
    ndf = ndf[np.abs(ndf) > thr].dropna()
    #print(ndf.shape)
    return ndf


def main():
    args = options()
    verbose = args.verbose
    biologID=pandas.read_csv(args.cpds, sep='\t')
    cpd2w=biologID.set_index('SEEDID')['Well'].T.to_dict()
    w2cpd=biologID.set_index('Well')['SEEDID'].T.to_dict()
    biologNames=pandas.read_csv(args.names, sep=',')
    w2name=biologNames.set_index('Well')['Compound'].T.to_dict()
    name2w=biologNames.set_index('Compound')['Well'].T.to_dict()
    w2class=biologNames.set_index('Well')['Class'].T.to_dict()

    keggID=pandas.read_csv(args.keggfile, sep='\t')
    seed2kegg=keggID.set_index('default')['KEGG'].T.to_dict()
    kegg2seed=keggID.set_index('KEGG')['default'].T.to_dict()

    names = ['Ammonia', 'L-Glutamic-Acid', 'L-Lysine', 'L-Serine', 'Urea']
    wells = ['A2', 'A12', 'B6', 'B10', 'A5']
    cpdsn = [w2cpd[w] for w in wells]

    fname='%s%s%s/%s' % (args.iopath, args.orglist, args.condition, args.fname)
    freq_thr = 0.0
    fluxthr = float(args.fluxthr)

    df = {}
    df_freq = {}
    df_norm = {}
    df_absnorm = {}
    for w in wells:
        cpdname=w2cpd[w]
        df[w]=pandas.read_csv(fname+'_fba_sol_'+cpdname+'.csv', index_col=0)
        #DataFrame.where(cond, other) Return an object of same shape as self and whose corresponding entries are from self where cond is True and otherwise are from other.
        df[w] = df[w].where(np.abs(df[w]) > fluxthr, 0)
        fdf = df[w].where(df[w] == 0, 1)
        fdf['freq'] = fdf[fdf==1].count(axis=1)/fdf.count(axis=1)
        df_freq[w] = fdf.copy()
        df_norm[w] = df[w].multiply(fdf['freq'], axis=0)
        df_absnorm[w] = np.abs(df_norm[w])
        df[w]['Mean'] = df[w].mean(axis=1, skipna=False)
        df_norm[w]['Mean'] = df_norm[w].mean(axis=1, skipna=False)
        df_absnorm[w]['Mean'] = df_absnorm[w].mean(axis=1, skipna=False)

    freq_df = mergeColumns(df_freq, 'freq', wells, freq_thr)
    means_df = mergeColumns(df, 'Mean', wells, 0)
    normmeans_df = mergeColumns(df_norm, 'Mean', wells, 0)
    absnormmeans_df = mergeColumns(df_absnorm, 'Mean', wells, 0)

    y = normmeans_df.divide(normmeans_df[wells[0]], axis=0)
    z = np.log2(np.abs(y))

    selthr = float(args.selthr)

    print('| %s | Same | Up | Down | ' % (args.orglist))
    print('| --- | --- | --- | --- |')
    for i in range(1,len(names)):
        cname = names[i][0:4]
        cwell = wells[i]
        vals = z[cwell]
        vals = vals.reset_index(level=0)
        vals['KEGG'] = vals['Row'].map(seed2kegg)
        #print(vals.shape)
        up = vals[vals[cwell] >= selthr].dropna()
        up.to_csv('%s%s%s/logfoldchange_amm_vs_%s_upFluxes.csv' % (args.iopath, args.orglist, args.condition, cname.lower()))
        #up['KEGG'].to_csv('%skegg_ids_%s%s_amm_vs_%s_up.csv' % (args.iopath, args.orglist, args.condition, cname.lower()), header=False, index=False)
        do = vals[vals[cwell] <= -1*selthr].dropna()
        do.to_csv('%s%s%s/logfoldchange_amm_vs_%s_doFluxes.csv' % (args.iopath, args.orglist, args.condition, cname.lower()))
        #do['KEGG'].to_csv('%skegg_ids_%s%s_amm_vs_%s_down.csv' % (args.iopath, args.orglist, args.condition, cname.lower()), header=False, index=False)
        sa = vals[np.abs(vals[cwell]) < selthr].dropna()
        sa.to_csv('%s%s%s/logfoldchange_amm_vs_%s_saFluxes.csv' % (args.iopath, args.orglist, args.condition, cname.lower()))
        #sa['KEGG'].to_csv('%skegg_ids_%s%s_amm_vs_%s_same.csv' % (args.iopath, args.orglist, args.condition, cname.lower()), header=False, index=False)
        print('| Ammonium vs %s | %d | %d | %d |' % (names[i], sa.shape[0], up.shape[0], do.shape[0]))
        formap = formatForKEGGMap(sa, up, do, '#7e7e7e %s #66C2A5' % (colors[names[i]]))
        formap.to_csv('%skegg_ids_%s%s_amm_vs_%s.csv' % (args.iopath, args.orglist, args.condition, cname.lower()), header=False, index=False, sep='\t')


def formatForKEGGMap(sa, up, do, colors):
    c = colors.split(' ')

    sa['Color'] = c[0]
    sa['Opacity'] = 0.8
    sa['Width'] = 'W15'
    df0 = sa[['KEGG', 'Color', 'Width', 'Opacity']]


    up['Color'] = c[1]
    up['Opacity'] = 0.8
    up['Width'] = 'W15'
    df1 = up[['KEGG', 'Color', 'Width', 'Opacity']]

    do['Color'] = c[2]
    do['Opacity'] = 0.8
    do['Width'] = 'W15'
    df2 = do[['KEGG', 'Color', 'Width', 'Opacity']]

    df = pandas.concat([df0, df1, df2])
    
    return df
    
def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-U', '--unmask', help='produce unmasked plots', action='store_true')
    parser.add_argument('-L', '--latex', help='print latex tab outputs', action='store_true')
    parser.add_argument('-M', '--markdown', help='print markdown tab outputs', action='store_true')
    parser.add_argument('-c', '--cpds', help='compunds file', default='../data/MPIRoots/singleNMedia/ncompounds.tsv')
    parser.add_argument('-k', '--keggfile', help='KEGG ID file', default='../data/ModelSEEDdata/KEGG.aliases.txt')
    parser.add_argument('-n', '--names', help='compunds wells name file', default='../data/MPIRoots/wellsNamesBiolog.csv')
    parser.add_argument('-f', '--fname', help='baseline file name', default='ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1')
    parser.add_argument('-p', '--iopath', help='path for input and output file', default='../outputs/')
    parser.add_argument('-o', '--orglist', help='list of organisms', default='Root9')
    parser.add_argument('-t', '--fluxthr', help='threshold below which flux is 0', default='1e-9')
    parser.add_argument('-s', '--selthr', help='threshold to select log fold change', default='1')
    parser.add_argument('-e', '--condition',  help='conditions e.g. excluding compounds',  default='_exclude_5N_and_nf_and_G12')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()

