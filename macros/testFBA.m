% Build a test ensemble
%
% Written by Matt Biggs, 2016

% This script has two purposes: to test the software environment and to
% serve as tutorial for using the software.
%
% TUTORIAL GUIDE:
% There are a few tasks that need to be accomplished before you can
% generate an ensemble and run a simulation: First, decide on a biomass
% objective function; Second, define media ingredients; Third, get a draft
% GENRE from the Model SEED for the organism and read it in; Fourth, set
% your parameters.
% After that, you're ready to reconstruct and simulate!

% Load universal reaction database
if ~exist('seed_rxns_mat','var')
    load 2018_seed_rxns
end

seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

GSMNMGenomicData = getGSMNMGenomeAnnotations('rhizobiumRoot491-reactions.tsv');

trimmed_rxnList = GSMNMGenomicData.rxn_GPR_mapping.rxns;
trimmed_gprs = GSMNMGenomicData.rxn_GPR_mapping.gprs;


biologicalData = getGSMNMGrowthConditions(seed_rxns_mat, 'growthMatrix_Root491.csv', 1);
%For test, subset the G/NG conditions
%Ammonium cpd00013
%Serine cpd00054
%growth_sources = [find(strcmp(seed_rxns_mat.mets, 'cpd00013')), find(strcmp(seed_rxns_mat.mets, 'cpd00054'))];
biologicalData.growthConditions = biologicalData.growthConditions(:,1:3);
biologicalData.nonGrowthConditions = biologicalData.nonGrowthConditions(:,1:2);

gs = biologicalData.growthXSources(1:3);
ngs = biologicalData.nonGrowthXSources(1:2);

% Store list of reaction indices to force inclusion during reconstruction
% process
Urxns2set = [find(ismember(seed_rxns_mat.rxns,trimmed_rxnList)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; % include spontaneous rxn05064
Uset2 = ones(size(Urxns2set));

% Include exchange reactions for all non-growth conditions, just so that
% it's the network itself--not the lack of exchange reactions--that prevents growth
% (n)gs == indexes of (non)growth sust sources
Xrxns2set = find(sum( abs(seed_rxns_mat.X([gs ngs],:)) ,1) > 0);
Xset2 = ones(size(Xrxns2set));

biologicalData.Urxns2set = Urxns2set;
biologicalData.Uset2 = Uset2;
biologicalData.Xrxns2set = Xrxns2set;
biologicalData.Xset2 = Xset2;


%------------------------------------------------------------------------
% Set parameters
%------------------------------------------------------------------------
params = struct;
params.sequential = 1;      % We want sequential gap filling (one growth condition at a time)
params.stochast = 1;        % Random element to gap filling
params.rndSeed = 1216;      % Setting the random seed ensures reproducibility
params.numModels2gen = 1;   % How many models do you want in your ensemble?
params.verbose = 1;         % 0 means we don't want to get updates about what's going on under the hood

%------------------------------------------------------------------------
% Gap fill a model!
%------------------------------------------------------------------------
tic
fprintf('Problem set up     ... success\n');
fprintf('Starting gap fill  ... (should finish in roughly 100 seconds)\n');
[modelList1] = build_network(seed_rxns_mat,biologicalData,params);
time2run = toc;
fprintf('Gap fill complete  ... success (%1.1f seconds)\n',time2run);
