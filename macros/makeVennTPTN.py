import argparse
import pandas
import numpy as np
import string
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import venn

def addMajorityCol(df):
    if 'Majority' in df.columns:
        print('Majority already in df, nothing to be done')
        return df
    tg = df[df==1].count(axis=1)
    tng = df[df==0].count(axis=1)
    df['TotG'] = tg
    df['TotNG'] = tng
    maj = df['TotG'] > df['TotNG']
    df['Majority'] = maj.where(~maj, 1)
    return

def getNetworkStats(g,ng):
    '''
    Compute accuracy, precision and recall for the network, excluding the conditions on which it was trained
    Perform column-wise operations on masked DF
    '''
    # 0 in g == False Negatives
    # 1 in g == True Positives
    # 0 in ng == True Negatives
    # 1 in ng == False Positives
    FN = g[g==0].count(axis=0)
    FP = ng[ng==1].count(axis=0)
    a = (TP+TN)/(TP+TN+FP+FN)
    p = (TP)/(TP+FP)
    r = (TP)/(TP+FN)
    return a,p,r


def main():
    args = options()
    verbose = args.verbose


    truesets = {}
    falsesets = {}

    for o in args.orglist.split(' '):
        tmpstr = args.fname.split('_')
            
        Nens = int(tmpstr[1])
        Ngcs = int(tmpstr[3])
        Nngcs = int(tmpstr[5])
        isStochW = int(tmpstr[8])

        fname=args.iopath+o+args.condition+'/'+args.fname
    
        dfc=pandas.read_csv(fname+'_conditions.csv')
        gdf=pandas.read_csv(fname+'_gc_tab.csv', index_col=0)
        ngdf=pandas.read_csv(fname+'_ngc_tab.csv', index_col=0)

        gcs = []
        ngcs = []
        for i in dfc.index:
            gcs.append(list(dfc.iloc[i,1:(1+Ngcs)].values))
            ngcs.append(list(dfc.iloc[i,(1+Ngcs):(1+Ngcs+Nngcs)].values))

        ngdf_mask = ngdf.copy()
        for i in range(len(ngcs)):
            ngdf_mask.iloc[:,i] = ~ngdf.iloc[:,i].index.isin(ngcs[i])
        ngdf_masked = ngdf.where(ngdf_mask, np.nan)

        gdf_mask = gdf.copy()
        for i in range(len(gcs)):
            gdf_mask.iloc[:,i] = ~gdf.iloc[:,i].index.isin(gcs[i])
        gdf_masked = gdf.where(gdf_mask, np.nan)

        gdf_masked_maj = gdf_masked.copy()
        addMajorityCol(gdf_masked_maj)
        ngdf_masked_maj = ngdf_masked.copy()
        addMajorityCol(ngdf_masked_maj)

        TP = set(gdf_masked_maj[gdf_masked_maj['Majority'] == 1].index)
        TN = set(ngdf_masked_maj[ngdf_masked_maj['Majority'] == 0].index)
        FP = set(gdf_masked_maj[gdf_masked_maj['Majority'] == 0].index)
        FN = set(ngdf_masked_maj[ngdf_masked_maj['Majority'] == 1].index)
        
        truesets[o] = TP ^ TN
        falsesets[o] = FP ^ FN


    alltptn = set.intersection(*list(truesets.values()))
    print('\n'.join(alltptn))
    venn.venn(truesets)
    plt.savefig(args.iopath+'venn_TN_TP.png')

    venn.venn(falsesets)
    plt.savefig(args.iopath+'venn_FN_FP.png')

    return truesets
    
def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-f', '--fname', help='baseline file name', default='ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1')
    parser.add_argument('-p', '--iopath', help='path for input and output file', default='../outputs/')
    parser.add_argument('-o', '--orglist', help='list of organisms', default='Root9 Root491 Root66D1')
    parser.add_argument('-c', '--condition',  help='conditions e.g. excluding compounds',  default='_exclude_not_found_and_G12')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
