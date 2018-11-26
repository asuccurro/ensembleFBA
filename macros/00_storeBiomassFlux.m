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

