%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/05/24     **
%**  last modified: 2018/06/04     **
%************************************
function [m] = runEnsemble(params)

% Load universal reaction database and add exchange rxns
load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

if params.ISTEST
  params.verbose = 1;
end

if ~exist('params.numGC','var')
    params.numGC = -1;
end

if ~exist('params.numNGC','var')
    params.numNGC = -1;
end

% Get the GSMNM data formatted with work with the SEED database
%       GSMNMData.biomassFn,growthXSources,growthConditions,nonGrowthXSources,nonGrowthConditions
[GSMNMData] = getGSMNMGrowthConditions(seed_rxns_mat, params.fileGrowthConditions, 1);

% Get the GSMNM gene-to-reaction mappings
%       GSMNMGenomicData.rxn_GPR_mapping
[GSMNMGenomicData] = getGSMNMGenomeAnnotations(params.fileAnnotations);
rxnList = GSMNMGenomicData.rxn_GPR_mapping.rxns;

% Force networks to contain reactions annotated from the genome
Urxns2set = [find(ismember(seed_rxns_mat.rxns,rxnList)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; % include spontaneous rxn05064
Uset2 = ones(size(Urxns2set));

% Include exchange reactions for all non-growth conditions, just so that
% it's the network itself--not the lack of exchange reactions--that prevents growth
Xrxns2set = find(sum( abs(seed_rxns_mat.X([GSMNMData.growthXSources(:); GSMNMData.nonGrowthXSources(:); GSMNMData.notForGapfillXSources(:)],:)) ,1) > 0);
Xset2 = ones(size(Xrxns2set));

GSMNMData.Urxns2set = Urxns2set;
GSMNMData.Uset2 = Uset2;
GSMNMData.Xrxns2set = Xrxns2set;
GSMNMData.Xset2 = Xset2;

full_growthConditions = GSMNMData.growthConditions;
full_nonGrowthConditions = GSMNMData.nonGrowthConditions;

ngc = params.numGC;
if ngc < 0
    ngc = floor(size(GSMNMData.growthConditions,2)/2);
end

nngc = params.numNGC;
if nngc < 0
    nngc = floor(size(GSMNMData.nonGrowthConditions,2)/2);
end

if params.ISTEST
  ntest=3;
  GSMNMData.growthConditions = GSMNMData.growthConditions(:,1:ntest);
  GSMNMData.nonGrowthConditions = GSMNMData.nonGrowthConditions(:,1:ntest);
  ngc = ntest;
  nngc = ntest;
  params.numModels2gen = ntest;
end
%------------------------------------------------------------------------
% Build a small ensemble!
%------------------------------------------------------------------------
GSMNMData.rxn_GPR_mapping = GSMNMGenomicData.rxn_GPR_mapping;
params.numGrowthConditions = ngc;
params.numNonGrowthConditions = nngc;
params.iterationThr = (ngc+nngc)*10;

fprintf('Starting build ensemble \n');
tic
[ensemble1] = build_ensemble(seed_rxns_mat,GSMNMData,params);
stseq1 = toc;

if length(ensemble1) == params.numModels2gen
    fprintf('Completed building ensemble  ... success\n');
else
    fprintf('Completed building ensemble  ... failure\n');
end

%------------------------------------------------------------------------
% Run eFBA!
%------------------------------------------------------------------------
tic
fprintf('Starting eFBA (should finish in roughly 1 second)\n');
[gc_growth] = ensembleFBA(ensemble1,seed_rxns_mat.Ex_names,full_growthConditions,0);
[ngc_growth] = ensembleFBA(ensemble1,seed_rxns_mat.Ex_names,full_nonGrowthConditions,0);
[xt_growth] = ensembleFBA(ensemble1,seed_rxns_mat.Ex_names,GSMNMData.notForGapfillConditions,0);
fprintf('eFBA run ... success\n');
stseq2 = toc;

m = struct;
m.ensemble = ensemble1;
m.gc_growth = gc_growth;
m.ngc_growth = ngc_growth;
m.xt_growth = xt_growth;
m.reconTime = stseq1;
m.solveTime = stseq2;
m.growthCpdList = seed_rxns_mat.mets(GSMNMData.growthXSources);
m.nonGrowthCpdList = seed_rxns_mat.mets(GSMNMData.nonGrowthXSources);
m.notForGapfillCpdList = seed_rxns_mat.mets(GSMNMData.notForGapfillXSources);

save(sprintf('%s/ensemble_%d_size_%d_gcs_%d_ngcs_stochasticWeights_%d.mat', params.fileOutPath, params.numModels2gen, ngc, nngc, params.stochast),'m');

end
