load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% load the ensembles of the 3 strains
mystrain = 'XXXORG';
mygenbank = 'XXXGID';
N=50;

load(fullfile('..', 'outputs', [mystrain '_exclude_not_found_and_G12'], 'ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1'))
[speciesGenomicData] = getGSMNMGenomeAnnotations(['MPI' mystrain '-reactions.tsv']);
ensemble = addGPRsGenBankIDs(m.ensemble,speciesGenomicData.rxn_GPR_mapping, mygenbank);

mediaCpd = 'XXXCPD';
blmedia = minimalmedia();
minimalMediaBase = zeros(length(seed_rxns_mat.mets),1);

x = [];
for k = 1:length(blmedia);
  x = [x, find(strcmp(seed_rxns_mat.mets, blmedia(k)))];
end

minimalMediaBase(x,1) = -100;
% Limiting nutrients
minimalMediaBase(find(strcmp(seed_rxns_mat.mets,  'cpd00009')), 1) = -5;
minimalMediaBase(find(strcmp(seed_rxns_mat.mets,  'cpd00027')), 1) = -5;
minimalMediaBase(find(strcmp(seed_rxns_mat.mets,  'cpd00048')), 1) = -5;

mediaXSources = [];
x = find(strcmp(seed_rxns_mat.mets, mediaCpd));
if length(x) > 0
    mediaXSources = [mediaXSources, x];
else
    fprintf(['AS-WARNING ' char(mediaCpd) ' (for media) not found in the rxn matrix\n']);
end

mediaConditions = repmat(minimalMediaBase,[1,length(mediaXSources)]);
for i = 1:length(mediaXSources)
    mediaConditions(mediaXSources(i),i) = -5;
end


curEnsemble = ensemble;
    
allGenes = cell(0,1);
for j = 1:N
    allGenes = [allGenes; curEnsemble{j}.genes];
end
allGenes = unique(allGenes);
    
geneEssentialityByNet = zeros(length(allGenes),N);

for l = 1:length(mediaXSources)
    for j = 1:N
        fprintf(['On network ', int2str(j), '\n'])
        curMod = curEnsemble{j};
        for k = 1:length(allGenes)
            curGene = allGenes{k};
            delMod = simulateGeneDeletion(curMod,curGene);
            delGrowth = fba_flex(delMod,seed_rxns_mat.Ex_names,mediaConditions(:,l),0);
            geneEssentialityByNet(k,j) = delGrowth < 1e-10;
        end
    end
    
    % Decide on essential genes
    essentialGeneIndicators = sum(geneEssentialityByNet,2) > N/2;
    essentialGenes = allGenes(essentialGeneIndicators > 0);
    %drugsHitEssentialGenes = genesWdrugMatches.drugIDs(ismember(genesWdrugMatches.genes,essentialGenes));
    
    fileName = ['../outputs/' mystrain '_exclude_not_found_and_G12/geneEssentiality_' mediaCpd '.mat']; 
    save(fileName,'geneEssentialityByNet','allGenes', 'essentialGenes');
end

N = sum(geneEssentialityByNet,2);
T=array2table(N,'RowNames',allGenes);
writetable(T,['../outputs/geneEssentiality/' mystrain '_XXXNAM.csv'],'WriteRowNames',true);
