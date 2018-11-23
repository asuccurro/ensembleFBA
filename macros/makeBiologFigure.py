
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
    TP = g[g==1].count(axis=0)
    TN = ng[ng==0].count(axis=0)
    FP = ng[ng==1].count(axis=0)
    a = (TP+TN)/(TP+TN+FP+FN)
    p = (TP)/(TP+FP)
    r = (TP)/(TP+FN)
    return a,p,r

def main():
    args = options()
    verbose = args.verbose
    biologID=pandas.read_csv(args.cpds, sep='\t')
    cpd2w=biologID.set_index('SEEDID')['Well'].T.to_dict()
    w2cpd=biologID.set_index('Well')['SEEDID'].T.to_dict()

    #ensemble_21_size_26_gcs_11_ngcs_stochasticWeights_0
    tmpstr = args.fname.split('_')

    Nens = int(tmpstr[1])
    Ngcs = int(tmpstr[3])
    Nngcs = int(tmpstr[5])
    isStochW = int(tmpstr[8])


    orgID = args.iopath.split('/')[-2].split('_')[0]
    ftit = '%s Growth vs No Growth predictions, Ensemble size %d' % (orgID, Nens)

    fname=args.iopath+args.fname
    dfc=pandas.read_csv(fname+'_conditions.csv')
    # Growth matrices on exp growth / no growth conditions
    # Columns = Networks
    # Rows = N sources
    #gdf=pandas.read_csv(fname+'_gc_growth.csv', header=None)
    #ngdf=pandas.read_csv(fname+'_ngc_growth.csv', header=None)
    gdf0=pandas.read_csv(fname+'_gc_tab.csv', index_col=0)
    ngdf=pandas.read_csv(fname+'_ngc_tab.csv', index_col=0)
    prdf=pandas.read_csv(fname+'_proteomics_growth.csv', index_col=0)

    # Add excluded conditions to DF
    # add prdf to gdf
    if 'exclude_A2-A5-A12-B6-B10' in fname or 'exclude_5N' in fname:
        gdf = pandas.concat([gdf0, prdf])
    else:
        gdf = gdf0
    
    #also generate random test dfs size=gdf.shape
    rndm_gdf = pandas.DataFrame(np.random.randint(0,2,size=(gdf.shape[0]-Ngcs,Nens)), columns=gdf.columns)
    rndm_ngdf = pandas.DataFrame(np.random.randint(0,2,size=(ngdf.shape[0]-Nngcs,Nens)), columns=ngdf.columns)

    # Matrix-like structure to link plate position with well label and compound SEED ID
    L=list(string.ascii_uppercase)
    pm_w = [12*[''] for n in range(8)]
    pm_cpd = [12*[''] for n in range(8)]
    for i in range(8):
        for j in range(12):
            pm_w[i][j] = '%s%d' % (L[i],j+1)
            pm_cpd[i][j] = w2cpd.get(pm_w[i][j],'NA')

    # ## Create masked DF to exclude the conditions used for trainings from stats computation
    ## Only gdf and ngdf, as prdf was not trained on conditions!
    gcs = []
    ngcs = []
    for i in dfc.index:
        gcs.append(list(dfc.iloc[i,1:(1+Ngcs)].values))
        ngcs.append(list(dfc.iloc[i,(1+Ngcs):(1+Ngcs+Nngcs)].values))

    ngdf_mask = ngdf.copy()
    for i in range(len(ngcs)):
        #print(ngdf.iloc[:,i].index.isin(ngcs[i]))
        ngdf_mask.iloc[:,i] = ~ngdf.iloc[:,i].index.isin(ngcs[i])
    ngdf_masked = ngdf.where(ngdf_mask, np.nan)

    gdf_mask = gdf.copy()
    for i in range(len(gcs)):
        gdf_mask.iloc[:,i] = ~gdf.iloc[:,i].index.isin(gcs[i])
    gdf_masked = gdf.where(gdf_mask, np.nan)

    #testng = ngdf_masked.iloc[:5,:3]
    #testg = gdf_masked.iloc[:8,:3]
    #addMajorityCol(testng)
    #addMajorityCol(testg)
    #ta,tp,tr = getNetworkStats(testg, testng)
    #print(ta['Majority'], tp['Majority'], tr['Majority'])

    gdf_masked_maj = gdf_masked.copy()
    addMajorityCol(gdf_masked_maj)
    ngdf_masked_maj = ngdf_masked.copy()
    addMajorityCol(ngdf_masked_maj)

    rndm_gdf_masked  = rndm_gdf.copy()
    rndm_ngdf_masked = rndm_ngdf.copy()
    rndm_gdf_masked  = rndm_gdf.where(gdf_mask, np.nan)
    rndm_ngdf_masked = rndm_ngdf.where(ngdf_mask, np.nan)
    addMajorityCol(rndm_gdf_masked)
    addMajorityCol(rndm_ngdf_masked)

    if args.unmask:
        plotBiologPlate(gdf, ngdf, args.unmask, fname+'_biologPlate_includingTrainingCond.png', ftit, pm_w, pm_cpd, prdf.index)
        addMajorityCol(gdf)
        addMajorityCol(ngdf)
        a,p,r = getNetworkStats(gdf,ngdf)
        addMajorityCol(rndm_gdf)
        addMajorityCol(rndm_ngdf)
        ra,rp,rr = getNetworkStats(rndm_gdf,rndm_ngdf)
        if args.latex:
            print('Unmasked Ensemble & %.3f & %.3f & %.3f \\\\' % (a['Majority'], p['Majority'], r['Majority']))
            print('Random Ensemble & %.3f & %.3f & %.3f \\\\' % (ra['Majority'], rp['Majority'], rr['Majority']))
        if args.markdown:
            print('| Unmasked Ensemble | %.3f | %.3f | %.3f | ' % (a['Majority'], p['Majority'], r['Majority']))
            print('| Random Ensemble | %.3f | %.3f | %.3f | ' % (ra['Majority'], rp['Majority'], rr['Majority']))
        return
    else:
        plotBiologPlate(gdf_masked_maj, ngdf_masked_maj, args.unmask, fname+'_biologPlate.png', ftit, pm_w, pm_cpd, prdf.index, args.markprot)
        a,p,r = getNetworkStats(gdf_masked_maj,ngdf_masked_maj)
        ra,rp,rr = getNetworkStats(rndm_gdf_masked,rndm_ngdf_masked)
        #addMajorityCol(rndm_gdf)
        #addMajorityCol(rndm_ngdf)
        #ra,rp,rr = getNetworkStats(rndm_gdf,rndm_ngdf)
        if args.latex:
            print('Masked Ensemble & %.3f & %.3f & %.3f \\\\' % (a['Majority'], p['Majority'], r['Majority']))
            print('Random Ensemble & %.3f & %.3f & %.3f \\\\' % (ra['Majority'], rp['Majority'], rr['Majority']))
        if args.markdown:
            print('| Masked Ensemble | %.3f | %.3f | %.3f | ' % (a['Majority'], p['Majority'], r['Majority']))
            print('| Random Ensemble | %.3f | %.3f | %.3f | ' % (ra['Majority'], rp['Majority'], rr['Majority']))
    return

def plotBiologPlate(gdf, ngdf, unmask, figname, figtit, pm_w, pm_cpd, pr_cpd, markprot=False):
    if unmask:
        if 'TotG' in list(gdf.columns):
            print('New columns TotG, TotNG and Majority should not be in the original DF!')
            return
    figPM = plt.figure(dpi=300)
    gsPM = gridspec.GridSpec(8, 12)
    colors_gng = ['#fa9fb5', '#e5f5f9']
    for i in range(8):
        for j in range(12):
            w = pm_w[i][j]
            cpd = pm_cpd[i][j]
            if cpd in list(gdf.index):
                ax = figPM.add_subplot(gsPM[i,j])
                #ax.get_xaxis().set_visible(False)
                #ax.get_yaxis().set_visible(False)
                if unmask:
                    ax.pie([gdf[gdf == 1].loc[cpd].count(), gdf[gdf==0].loc[cpd].count()], colors=colors_gng, startangle=90)
                else:
                    ax.pie([gdf.loc[cpd,'TotG'], gdf.loc[cpd,'TotNG']], colors=colors_gng, startangle=90)
                ax.pie([1, 0], colors=colors_gng,radius=0.5,startangle=90)
                ax.text(-1, 1.25, w, fontsize=6)
                #ax.set_facecolor('#fa9fb5')
                #gdf_masked.loc[cpd].plot.hist()
            elif cpd in list(ngdf.index):
                ax = figPM.add_subplot(gsPM[i,j])
                #ax.get_xaxis().set_visible(False)
                #ax.get_yaxis().set_visible(False)
                if unmask:
                    ax.pie([ngdf[ngdf == 1].loc[cpd].count(), ngdf[ngdf==0].loc[cpd].count()], colors=colors_gng, startangle=90)
                else:
                    ax.pie([ngdf.loc[cpd,'TotG'], ngdf.loc[cpd,'TotNG']], colors=colors_gng, startangle=90)
                ax.pie([0, 1], colors=colors_gng,radius=0.5,startangle=90)
                ax.text(-1, 1.25, w, fontsize=6)
            else:
                ax = figPM.add_subplot(gsPM[i,j])
                #ax.axis('off')
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.pie([0, 1], colors=['w', 'w'],radius=0.5,startangle=90)
                ax.text(-1, 1.25, w, fontsize=6)
                if i==0 and j==0:
                    ax.text(0., 0., 'control', horizontalalignment='center', fontsize=6)
                else:
                    ax.text(0., 0., 'N/A', horizontalalignment='center', fontsize=6)
            if cpd in pr_cpd and markprot:
                autoAxis = ax.axis()
                rec = matplotlib.patches.Rectangle((autoAxis[0],autoAxis[2]+0.1),(autoAxis[1]-autoAxis[0]),(autoAxis[3]-autoAxis[2])+0.5,fill=False,lw=1, color='red')
                rec = ax.add_patch(rec)
                rec.set_clip_on(False)
    figPM.suptitle(figtit)
    figPM.savefig(figname)

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-P', '--markprot', help='mark the proteomics N sources', action='store_true')
    parser.add_argument('-U', '--unmask', help='produce unmasked plots', action='store_true')
    parser.add_argument('-L', '--latex', help='print latex tab outputs', action='store_true')
    parser.add_argument('-M', '--markdown', help='print markdown tab outputs', action='store_true')
    parser.add_argument('-c', '--cpds', help='compunds file', default='../data/MPIRoots/singleNMedia/ncompounds.tsv')
    parser.add_argument('-f', '--fname', help='baseline file name', default='ensemble_21_size_26_gcs_11_ngcs_stochasticWeights_0')
    parser.add_argument('-p', '--iopath', help='path for input and output file', default='../outputs/')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
