%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
%** Template to extract the biomass fluxes
%************************************
ensembleFname='XXXFNAME';
outpath='../outputs/XXXDNAME/';
% Load the struct containing the ensemble
load(fullfile('..', 'outputs', 'XXXDNAME', ensembleFname))
N=[m.gc_growth; m.ngc_growth];
TG=array2table(N,'RowNames', [m.growthCpdList; m.nonGrowthCpdList]);
writetable(TG,[outpath ensembleFname, '_biomassFluxes.csv'], 'WriteRowNames',true);

load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

allCpds=[m.growthCpdList; m.nonGrowthCpdList; m.notForGapfillCpdList];
biomass_fluxes = computeFBAsol(m.ensemble, seed_rxns_mat, allCpds, ensembleFname, outpath, 0);
%N=biomass_fluxes;
%TG=array2table(N,'RowNames',allCpds);
%writetable(TG,[outpath ensembleFname, '_fba_growth.csv'], 'WriteRowNames',true);
