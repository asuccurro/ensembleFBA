%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
%** Analyse the ensemble for the Root491 draft GSMN
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

%dlmwrite([outpath ensembleFname, '_gc_growth.csv'], m.gc_growth>0, ',');
%dlmwrite([outpath ensembleFname, '_ngc_growth.csv'], m.ngc_growth>0, ',');

N=uint8(m.xt_growth>0);
TG=array2table(N,'RowNames',m.notForGapfillCpdList);
writetable(TG,[outpath ensembleFname, '_proteomics_growth.csv'], 'WriteRowNames',true);

%dlmwrite([outpath ensembleFname, '_cond.csv'], c, ',');

e = m.ensemble;
r = seed_rxns_mat.Ex_names;
c = m.ensemble{1}.notForGapfillConditions;

% for every media condition and for every network store the solution fluxes
% store in place 1 the reactions, in place 2 the EX_rxns
solutions = cell(size(c,2),2);
for i = 1:size(c,2)
  s_rx = zeros(size(seed_rxns_mat.rxns,1),length(e));
  s_ex = zeros(size(seed_rxns_mat.Ex_names,1),length(e));
  for j = 1:length(e)
    % match model reactions to matrix; distinguish exchange from rxns as are named differently
    model = e{j};
    rxnindex = [];
    exrindex = [];
    rr = [];
    ee = [];
    for z = 1:length(model.rxns)
      ef = find(strcmp(seed_rxns_mat.Ex_names, model.rxns(z)));
      rf = find(strcmp(seed_rxns_mat.rxns, model.rxns(z)));
      if length(rf) > 0
        rxnindex = [rxnindex rf];
        rr = [rr z];
      elseif length(ef)>0
        exrindex = [exrindex ef];
        ee = [ee z];
      end
    end
    [growth,x] = fba_flex(model,r,c(:,i),1);
    s_rx(rxnindex,j) = x(rr);
    s_ex(exrindex,j) = x(ee);
  end
  solutions{i,1} = s_rx;
  solutions{i,2} = s_ex;
  % Actually we do not care about Ex reactions
  N = s_rx;
  T=array2table(N,'RowNames',seed_rxns_mat.rxns);
  writetable(T,[outpath ensembleFname '_fba_sol_' char(m.notForGapfillCpdList(i)) '.csv'], 'WriteRowNames',true);
  dlmwrite([outpath ensembleFname, '_exc_rxns_' char(m.notForGapfillCpdList(i)) '.csv'], s_ex, ',');
end

% To do in python
  % 1. Check "enrichment" of reactions in ensemble
  % 2. Average fluxes and make heat map



%
