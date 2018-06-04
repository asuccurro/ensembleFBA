%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
%** Analyse the ensemble for the Root491 draft GSMN
%************************************
ensembleFname='Ensemble_5_gcs_5_ngcs';
outpath='../outputs/ensemble/';
% Load universal reaction database and add exchange rxns
load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Load the struct containing the ensemble
load(ensembleFname)
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

dlmwrite([outpath ensembleFname, '_gc_growth.csv'], m.gc_growth>0, ',');
dlmwrite([outpath ensembleFname, '_ngc_growth.csv'], m.ngc_growth>0, ',');
