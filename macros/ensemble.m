%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/05/24     **
%**  last modified: 2018/05/24     **
%************************************
%** Gapfill sequentially the Root491 draft GSMN
%** Based on the scripts from Matt Biggs, 2016
%************************************

% Load universal reaction database and add exchange rxns
load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% Get the GSMNM data formatted with work with the SEED database
%       GSMNMData.biomassFn,growthXSources,growthConditions,nonGrowthXSources,nonGrowthConditions
[GSMNMData] = getGSMNMGrowthConditions(seed_rxns_mat, 'growthMatrix_Root491.csv', 1);

% Get the GSMNM gene-to-reaction mappings
%       GSMNMGenomicData.rxn_GPR_mapping
[GSMNMGenomicData] = getGSMNMGenomeAnnotations('rhizobiumRoot491-reactions.tsv');
rxnList = GSMNMGenomicData.rxn_GPR_mapping.rxns;

% Force networks to contain reactions annotated from the genome
Urxns2set = [find(ismember(seed_rxns_mat.rxns,rxnList)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; % include spontaneous rxn05064
Uset2 = ones(size(Urxns2set));

% Include exchange reactions for all non-growth conditions, just so that
% it's the network itself--not the lack of exchange reactions--that prevents growth
Xrxns2set = find(sum( abs(seed_rxns_mat.X([GSMNMData.growthXSources(:); GSMNMData.nonGrowthXSources(:)],:)) ,1) > 0);
Xset2 = ones(size(Xrxns2set));

% Set parameters
GSMNMData.Urxns2set = Urxns2set;
GSMNMData.Uset2 = Uset2;
GSMNMData.Xrxns2set = Xrxns2set;
GSMNMData.Xset2 = Xset2;

params = struct;
params.sequential = 1;
params.stochast = 0;
params.numModels2gen = 1;
params.verbose = 1;

full_growthConditions = GSMNMData.growthConditions;
full_nonGrowthConditions = GSMNMData.nonGrowthConditions;
GSMNMData.growthConditions = GSMNMData.growthConditions(:,1:3);
GSMNMData.nonGrowthConditions = GSMNMData.nonGrowthConditions(:,1:3);
%------------------------------------------------------------------------
% Build a small ensemble!
%------------------------------------------------------------------------
fprintf('Starting build ensemble    (should finish in roughly 50 seconds)\n');
GSMNMData.rxn_GPR_mapping = GSMNMGenomicData.rxn_GPR_mapping;
params.fractionUrxns2set = 0.8;
params.rndSequence = 1;
params.numModels2gen = 3;
params.numGrowthConditions = 3; %floor(size(GSMNMData.growthConditions,2)/2);
params.numNonGrowthConditions = 3; %floor(size(GSMNMData.nonGrowthConditions,2)/2);
[ensemble1] = build_ensemble(seed_rxns_mat,GSMNMData,params);

if length(ensemble1) == 3
    fprintf('Completed building ensemble  ... success\n');
else
    fprintf('Completed building ensemble  ... failure\n');
end

%------------------------------------------------------------------------
% Run eFBA!
%------------------------------------------------------------------------
fprintf('Starting eFBA (should finish in roughly 1 second)\n');
[gc_growth] = ensembleFBA(ensemble1,seed_rxns_mat.Ex_names,full_growthConditions,0);
[ngc_growth] = ensembleFBA(ensemble1,seed_rxns_mat.Ex_names,full_nonGrowthConditions,0);
fprintf('eFBA run ... success\n');
