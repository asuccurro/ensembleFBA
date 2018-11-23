import pandas
import numpy as np

fname='../outputs/Root66D1_exclude_5N_and_nf_and_G12/ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1'

tmpstr = fname.split('/')[-1].split('_')

Nens = int(tmpstr[1])
Ngcs = int(tmpstr[3])
Nngcs = int(tmpstr[5])
isStochW = int(tmpstr[8])


dfc=pandas.read_csv(fname+'_conditions.csv')
gdf0=pandas.read_csv(fname+'_gc_tab.csv', index_col=0)
ngdf=pandas.read_csv(fname+'_ngc_tab.csv', index_col=0)
prdf=pandas.read_csv(fname+'_proteomics_growth.csv', index_col=0)

if 'exclude_A2-A5-A12-B6-B10' in fname or 'exclude_5N' in fname:
    gdf = pandas.concat([gdf0, prdf])
else:
    gdf = gdf0



rndm_gdf = pandas.DataFrame(np.random.randint(0,2,size=(gdf.shape[0]-Ngcs,Nens)), columns=gdf.columns)
rndm_ngdf = pandas.DataFrame(np.random.randint(0,2,size=(ngdf.shape[0]-Nngcs,Nens)), columns=ngdf.columns)

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
    
rndm_gdf_masked  = rndm_gdf.copy()
rndm_ngdf_masked = rndm_ngdf.copy()
rndm_gdf_masked  = rndm_gdf.where(gdf_mask, np.nan)
rndm_ngdf_masked = rndm_ngdf.where(ngdf_mask, np.nan)
addMajorityCol(rndm_gdf_masked)
addMajorityCol(rndm_ngdf_masked)
