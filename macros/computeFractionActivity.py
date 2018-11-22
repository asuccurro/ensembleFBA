#!/usr/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/24     **
#**  last modified: 2018/05/24     **
#************************************

import argparse
import pandas
import numpy as np
import string
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from makeBiologFigure import getNetworkStats, addMajorityCol

def main():
    args = options()
    verbose = args.verbose
    biologID=pandas.read_csv(args.cpds, sep='\t')
    cpd2w=biologID.set_index('SEEDID')['Well'].T.to_dict()
    w2cpd=biologID.set_index('Well')['SEEDID'].T.to_dict()
    biologNames=pandas.read_csv(args.names, sep=',')
    w2name=biologNames.set_index('Well')['Compound'].T.to_dict()
    w2class=biologNames.set_index('Well')['Class'].T.to_dict()

    print('\n****N.B these number do not take into account excluded conditions!***\n')
    
    #ensemble_21_size_26_gcs_11_ngcs_stochasticWeights_0
    tmpstr = args.fname.split('_')

    Nens = int(tmpstr[1])
    Ngcs = int(tmpstr[3])
    Nngcs = int(tmpstr[5])
    isStochW = int(tmpstr[8])

    #Root9 Root491 Root66D1
    orglist = args.orglist.split(' ')

    flist = {}
    conditions_df = {}
    growth_df = {}
    nongrowth_df = {}
    ngdf_masked = {}
    gdf_masked = {}
    for o in orglist:
        flist[o] = args.iopath+o+args.condition+'/'+args.fname
        conditions_df[o] = pandas.read_csv(flist[o]+'_conditions.csv')
        growth_df[o] = pandas.read_csv(flist[o]+'_gc_tab.csv', index_col=0)
        nongrowth_df[o] = pandas.read_csv(flist[o]+'_ngc_tab.csv', index_col=0)
        
        gcs = []
        ngcs = []
        for i in conditions_df[o].index:
            gcs.append(list(conditions_df[o].iloc[i,1:(1+Ngcs)].values))
            ngcs.append(list(conditions_df[o].iloc[i,(1+Ngcs):(1+Ngcs+Nngcs)].values))

        ngdf_mask = nongrowth_df[o].copy()
        for i in range(len(ngcs)):
            #print(ngdf.iloc[:,i].index.isin(ngcs[i]))
            ngdf_mask.iloc[:,i] = ~nongrowth_df[o].iloc[:,i].index.isin(ngcs[i])
        ngdf_masked[o] = nongrowth_df[o].where(ngdf_mask, np.nan)

        gdf_mask = growth_df[o].copy()
        for i in range(len(gcs)):
            gdf_mask.iloc[:,i] = ~growth_df[o].iloc[:,i].index.isin(gcs[i])
        gdf_masked[o] = growth_df[o].where(gdf_mask, np.nan)
        addMajorityCol(gdf_masked[o])
        addMajorityCol(ngdf_masked[o])

    with open(args.iopath+'activity'+args.condition+'.csv', 'w') as ofile:
        ofile.write('Well,Name,Class,%s\n' % ','.join(orglist))
        for c in list(gdf_masked[orglist[0]].index)+list(ngdf_masked[orglist[0]].index):
            cstr = '%s,%s,%s' % (cpd2w[c], w2name[cpd2w[c]], w2class[cpd2w[c]])
            for o in orglist:
                if c in gdf_masked[o].index:
                    activity = gdf_masked[o].loc[c,'TotG']/(gdf_masked[o].loc[c,'TotG'] + gdf_masked[o].loc[c,'TotNG'])
                elif c in ngdf_masked[o].index:
                    activity = ngdf_masked[o].loc[c,'TotG']/(ngdf_masked[o].loc[c,'TotG'] + ngdf_masked[o].loc[c,'TotNG'])
                else:
                    print('AAAA')
                cstr = cstr+ (',%.3f' % activity)
            ofile.write('%s\n' % cstr)
    

    return gdf_masked,ngdf_masked



def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-U', '--unmask', help='produce unmasked plots', action='store_true')
    parser.add_argument('-L', '--latex', help='print latex tab outputs', action='store_true')
    parser.add_argument('-M', '--markdown', help='print markdown tab outputs', action='store_true')
    parser.add_argument('-c', '--cpds', help='compunds file', default='../data/MPIRoots/singleNMedia/ncompounds.tsv')
    parser.add_argument('-n', '--names', help='compunds wells name file', default='../data/MPIRoots/wellsNamesBiolog.csv')
    parser.add_argument('-f', '--fname', help='baseline file name', default='ensemble_50_size_21_gcs_10_ngcs_stochasticWeights_1')
    parser.add_argument('-p', '--iopath', help='path for input and output file', default='../outputs/')
    parser.add_argument('-o', '--orglist', help='list of organisms', default='Root9 Root491 Root66D1')
    parser.add_argument('-e', '--condition',  help='conditions e.g. excluding compounds',  default='_exclude_not_found_and_G12')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
