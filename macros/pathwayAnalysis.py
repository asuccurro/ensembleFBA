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


def mergeColumns(dfdic, coln, cols, thr=0):
    ndf = pandas.DataFrame(dfdic[cols[0]][coln])
    ndf.columns = [cols[0]]
    for c in cols[1:]:
        ndf[c] = dfdic[c][coln]
    print(ndf.shape)
    ndf = ndf[np.abs(ndf) > thr].dropna()
    print(ndf.shape)
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

    growth_df = pandas.read_csv(fname+'_fba_growth.csv', index_col=0)

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
    
    for i in range(1,len(names)):
        cname = names[i][0:4]
        cwell = wells[i]
        vals = z[cwell]
        vals = vals.reset_index(level=0)
        vals['KEGG'] = vals['Row'].map(seed2kegg)
        print(vals.head())
        up = vals.dropna()[vals[cwell] >= selthr]
        up.to_csv('%s%s%s/logfoldchange_upFluxes.csv' % (args.iopath, args.orglist, args.condition))
        do = vals.dropna()[vals[cwell] <= -1*selthr]
        do.to_csv('%s%s%s/logfoldchange_doFluxes.csv' % (args.iopath, args.orglist, args.condition))
        sa = vals.dropna()[np.abs(vals[cwell]) < selthr]
        sa.to_csv('%s%s%s/logfoldchange_saFluxes.csv' % (args.iopath, args.orglist, args.condition))
        


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

