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
    growth_fl = {}
    gf_masked = {}
    gng_masked = {}
    gf_means = {}
    for o in orglist:
        flist[o] = args.iopath+o+args.condition+'/'+args.fname
        conditions_df[o] = pandas.read_csv(flist[o]+'_conditions.csv')
        growth_fl[o] = pandas.read_csv(flist[o]+'_biomassFluxes.csv', index_col=0)
        growth_fl[o] = growth_fl[o].where(growth_fl[o] > 0.000000001, 0)
        
        gfc = []
        for i in conditions_df[o].index:
            gfc.append(list(conditions_df[o].iloc[i,1:(1+Ngcs+Nngcs)].values))

        gf_mask = growth_fl[o].copy()
        for i in range(len(gfc)):
            #print(ngdf.iloc[:,i].index.isin(ngcs[i]))
            gf_mask.iloc[:,i] = ~growth_fl[o].iloc[:,i].index.isin(gfc[i])
            
        gf_masked[o] = growth_fl[o].where(gf_mask, np.nan)
        gf_means[o] = gf_masked[o].mean(1)
        
        gng_masked[o] = growth_fl[o].where(growth_fl[o] < 0.000000001, 1)
        gng_masked[o] = gng_masked[o].where(gf_mask, np.nan)

        addMajorityCol(gng_masked[o])

    with open(args.iopath+'activity_as_growth_fraction'+args.condition+'.csv', 'w') as ofile:
        ofile.write('Well,Name,Class,%s\n' % ','.join(orglist))
        for c in list(gng_masked[orglist[0]].index):
            cstr = '%s,%s,%s' % (cpd2w[c], w2name[cpd2w[c]], w2class[cpd2w[c]])
            for o in orglist:
                if c in gng_masked[o].index:
                    activity = gng_masked[o].loc[c,'TotG']/(gng_masked[o].loc[c,'TotG'] + gng_masked[o].loc[c,'TotNG'])
                else:
                    print('AAAA')
                cstr = cstr+ (',%.3f' % activity)
            ofile.write('%s\n' % cstr)

    with open(args.iopath+'activity_as_growth_average'+args.condition+'.csv', 'w') as ofile:
        ofile.write('Well,Name,Class,%s\n' % ','.join(orglist))
        for c in list(gf_means[orglist[0]].index):
            cstr = '%s,%s,%s' % (cpd2w[c], w2name[cpd2w[c]], w2class[cpd2w[c]])
            for o in orglist:
                if c in gf_means[o].index:
                    activity = gf_means[o][c]
                else:
                    print('AAAA')
                cstr = cstr+ (',%.3f' % activity)
            ofile.write('%s\n' % cstr)

    with open(args.iopath+'activity_as_weighted_growth_average'+args.condition+'.csv', 'w') as ofile:
        ofile.write('Well,Name,Class,%s\n' % ','.join(orglist))
        for c in list(gf_means[orglist[0]].index):
            cstr = '%s,%s,%s' % (cpd2w[c], w2name[cpd2w[c]], w2class[cpd2w[c]])
            for o in orglist:
                if c in gf_means[o].index:
                    activity = gf_means[o][c]*(gng_masked[o].loc[c,'TotG']/(gng_masked[o].loc[c,'TotG'] + gng_masked[o].loc[c,'TotNG']))
                else:
                    print('AAAA')
                cstr = cstr+ (',%.3f' % activity)
            ofile.write('%s\n' % cstr)

    return gf_masked



def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-U', '--unmask', help='produce unmasked plots', action='store_true')
    parser.add_argument('-L', '--latex', help='print latex tab outputs', action='store_true')
    parser.add_argument('-M', '--markdown', help='print markdown tab outputs', action='store_true')
    parser.add_argument('-c', '--cpds', help='compunds file', default='../data/MPIRoots/singleNMedia/ncompounds.tsv')
    parser.add_argument('-n', '--names', help='compunds wells name file', default='../data/MPIRoots/wellsNamesBiolog.csv')
    parser.add_argument('-f', '--fname', help='baseline file name', default='ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1')
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
