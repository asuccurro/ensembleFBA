%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
%** Template to analyse ensemble
%************************************
ensembleFname='XXXFNAME';
outpath='../outputs/XXXDNAME/';
% Load universal reaction database and add exchange rxns
load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Load the struct containing the ensemble
load(fullfile('..', 'outputs', 'XXXDNAME', ensembleFname))
% m =
%
%   struct with fields:
%
%       ensemble: {21×1 cell}
%      gc_growth: [53×21 double]
%     ngc_growth: [22×21 double]
%      reconTime: 19.8216
%      solveTime: 70.4936

% Check the growth / non growth conditions used, store in a table
fid=fopen([outpath ensembleFname, '_conditions.csv'],'w');
header = ['Network,'];
for i=1:size(m.ensemble{1}.growthConditions,2)
  header = [header sprintf('GC%d,', i)];
end
for i=1:size(m.ensemble{1}.nonGrowthConditions,2)
  header = [header sprintf('NGC%d,', i)];
end
fprintf(fid, [header 'AS\n']);
for i=1:length(m.ensemble)
  reconCon = getReconConditions(seed_rxns_mat, m.ensemble{i});
  fprintf(fid, [num2str(i) ',' sprintf('%s,', reconCon.growth{:}) sprintf('%s,', reconCon.nonGrowth{:}) 'AS\n']);
end
fclose(fid);

N=uint8(m.gc_growth>0);
TG=array2table(N,'RowNames',m.growthCpdList);
writetable(TG,[outpath ensembleFname, '_gc_tab.csv'], 'WriteRowNames',true);

N=uint8(m.ngc_growth>0);
TNG=array2table(N,'RowNames',m.nonGrowthCpdList);
writetable(TNG,[outpath ensembleFname, '_ngc_tab.csv'], 'WriteRowNames',true);

N=uint8(m.xt_growth>0);
TG=array2table(N,'RowNames',m.notForGapfillCpdList);
writetable(TG,[outpath ensembleFname, '_nfg_tab.csv'], 'WriteRowNames',true);

N=[m.gc_growth; m.ngc_growth; m.xt_growth];
TG=array2table(N,'RowNames', [m.growthCpdList; m.nonGrowthCpdList; m.notForGapfillCpdList]);
writetable(TG,[outpath ensembleFname, '_biomass_tab.csv'], 'WriteRowNames',true);

allCpds=[m.growthCpdList; m.nonGrowthCpdList; m.notForGapfillCpdList];
biomass_fluxes = computeFBAsol(m.ensemble, seed_rxns_mat, allCpds, 'allCond', ensembleFname, outpath, 0, 0);


