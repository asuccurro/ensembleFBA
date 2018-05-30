

load seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

pa = getPA14GrowthConditions(seed_rxns_mat);
pas = pa.biomassFn < 0;
pap = pa.biomassFn > 0;
bm_substrates = seed_rxns_mat.mets(pas);
bm_products = seed_rxns_mat.mets(pap);

printGrepLoop(bm_substrates)
printGrepLoop(bm_products)

lb_minimal_base = zeros(length(seed_rxns_mat.mets),1);
lbm = [4926,2213,3321,6862,7930,6753,6714,2938,7928,7934,7222];
lb_minimal_base(lbm,1) = -1000;

minimalmedia = seed_rxns_mat.mets(lbm);

printGrepLoop(minimalmedia)

growth_carbon_sources = [7578,1345,6524,980];
nongrowth_carbon_sources = [8832,1859,8289];

gcs = seed_rxns_mat.mets(growth_carbon_sources);
ngcs = seed_rxns_mat.mets(nongrowth_carbon_sources);

printGrepLoop(gcs)
printGrepLoop(ngcs)

n = length(growth_carbon_sources);
growthConditions = repmat(lb_minimal_base,[1,n]);
for i = 1:n
    growthConditions(growth_carbon_sources(i),i) = -10;
end

nonGrowthConditions = repmat(lb_minimal_base,[1,length(nongrowth_carbon_sources)]);
for i = 1:length(nongrowth_carbon_sources)
    nonGrowthConditions(nongrowth_carbon_sources(i),i) = -10;
end

%%%%%%%%% end of Matt's way

%%%%%%%%% start my way, using old seed_rxns_mat
biomassFn = myBiomassBuilder(seed_rxns_mat, bm_substrates, bm_products);
isequal(biomassFn, pa.biomassFn)

[minimalMediaBase, growthConditions2, nonGrowthConditions2] = myMediaBuilder(seed_rxns_mat, minimalmedia, gcs, ngcs);

isequal(lb_minimal_base, minimalMediaBase)
isequal(growthConditions, growthConditions2)
isequal(nonGrowthConditions, nonGrowthConditions2)

%%%%%%%%% until here everything checks


%%%%%%%%% Matt's code for model recon

%%%%%%%%% GPR map, independent of DB
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

Urxns2set = [find(ismember(seed_rxns_mat.rxns,trimmed_rxnList)); find(ismember(seed_rxns_mat.rxns,'rxn05064'))]; 
Uset2 = ones(size(Urxns2set));

Xrxns2set = find(sum( abs(seed_rxns_mat.X([growth_carbon_sources nongrowth_carbon_sources],:)) ,1) > 0);
Xset2 = ones(size(Xrxns2set));

biologicalData = struct; % Store data as Matlab struct 
biologicalData.growthConditions = growthConditions;
biologicalData.nonGrowthConditions = nonGrowthConditions;
biologicalData.biomassFn = pa.biomassFn;
biologicalData.Urxns2set = Urxns2set;
biologicalData.Uset2 = Uset2;
biologicalData.Xrxns2set = Xrxns2set;
biologicalData.Xset2 = Xset2;



%%%%%%%%%% Move to new DB
old_seed_rxns_mat = seed_rxns_mat;


load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

new_seed_rxns_mat = seed_rxns_mat;

biomassFn = myBiomassBuilder(seed_rxns_mat, bm_substrates, bm_products);

biomassFn_PA14 = biomassFn;

save('biomassFn_PA14_2018.mat','biomassFn_PA14');

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
params.verbose = 2;         % 0 means we don't want to get updates about what's going on under the hood

tic
fprintf('Problem set up     ... success\n');
fprintf('Starting gap fill  ... (should finish in roughly 100 seconds)\n');
[modelList1] = build_network(old_seed_rxns_mat,biologicalData,params);
time2run = toc;
fprintf('Gap fill complete  ... success (%1.1f seconds)\n',time2run);

tic
fprintf('Problem set up     ... success\n');
fprintf('Starting gap fill  ... (should finish in roughly 100 seconds)\n');
[modelList1] = build_network(new_seed_rxns_mat,myData,params);
time2run = toc;
fprintf('Gap fill complete  ... success (%1.1f seconds)\n',time2run);
