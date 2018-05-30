load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

pa = getPA14GrowthConditions(seed_rxns_mat);
pas = pa.biomassFn < 0;
pap = pa.biomassFn > 0;
bm_substrates = seed_rxns_mat.mets(pas);
bm_products = seed_rxns_mat.mets(pap);

lb_minimal_base = zeros(length(seed_rxns_mat.mets),1);
lbm = [4926,2213,3321,6862,7930,6753,6714,2938,7928,7934,7222];
lb_minimal_base(lbm,1) = -1000;

minimalmedia = seed_rxns_mat.mets(lbm);

growth_carbon_sources = [7578,1345,6524,980];
nongrowth_carbon_sources = [8832,1859,8289];

gcs = seed_rxns_mat.mets(growth_carbon_sources);
ngcs = seed_rxns_mat.mets(nongrowth_carbon_sources);

%%%%%%%%% Matt's code for model recon
fid = fopen('PA14_reference_genome.rxntbl','r');
rxns = textscan(fid, '%s%s%s%s%s','Delimiter','\t');
fclose(fid);

trimmed_rxnList = char(rxns{1});
trimmed_rxnList = trimmed_rxnList(:,1:8);
gprs = strfind(rxns{5},'fig');
rxns_with_genomic_evidence = ~cellfun(@isempty,gprs);
trimmed_rxnList = cellstr(trimmed_rxnList(rxns_with_genomic_evidence,:));
trimmed_gprs = rxns{5};
trimmed_gprs = cellfun(@(orig,old,new) strrep(orig,old,new), trimmed_gprs, repmat({'fig|208963.12.'},[length(trimmed_gprs),1]), repmat({''},[length(trimmed_gprs),1]),'UniformOutput',false);
trimmed_gprs = trimmed_gprs(rxns_with_genomic_evidence);

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

%%%%%%%%%% Move to new DB
old_seed_rxns_mat = seed_rxns_mat;

load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

new_seed_rxns_mat = seed_rxns_mat;

biomassFn = myBiomassBuilder(seed_rxns_mat, bm_substrates, bm_products);

[minimalMediaBase, growthConditions2, nonGrowthConditions2, gid, ngid] = myMediaBuilder(seed_rxns_mat, minimalmedia, gcs, ngcs);

Urxns2newset = [find(ismember(seed_rxns_mat.rxns,trimmed_rxnList)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; 
Unewset2 = ones(size(Urxns2newset));

Xrxns2newset = find(sum( abs(seed_rxns_mat.X([gid ngid],:)) ,1) > 0);
Xnewset2 = ones(size(Xrxns2newset));

myData = struct; % Store data as Matlab struct 
myData.growthConditions = growthConditions2;
myData.nonGrowthConditions = nonGrowthConditions2;
myData.biomassFn = biomassFn;
myData.Urxns2set = Urxns2newset;
myData.Uset2 = Unewset2;
myData.Xrxns2set = Xrxns2newset;
myData.Xset2 = Xnewset2;

params = struct;
params.sequential = 1;      % We want sequential gap filling (one growth condition at a time)
params.stochast = 1;        % Random element to gap filling
params.rndSeed = 1216;      % Setting the random seed ensures reproducibility
params.numModels2gen = 1;   % How many models do you want in your ensemble?
params.verbose = 1;         % 0 means we don't want to get updates about what's going on under the hood

tic
fprintf('Problem set up     ... success\n');
fprintf('Starting gap fill  ... (should finish in roughly 100 seconds)\n');
[modelList1] = build_network(new_seed_rxns_mat,myData,params);
time2run = toc;
fprintf('Gap fill complete  ... success (%1.1f seconds)\n',time2run);
