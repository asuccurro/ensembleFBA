
load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

% load the ensembles of the 3 strains
ensembles = cell(50,3);
speciesOrder = {'Root9', 'Root491', 'Root66D1'};
speciesFig = {'fig|1736604.3.', 'fig|1736548.3.', 'fig|1736582.3.'};
speciesGenBank = {'ASE33_', 'ASD46_', 'ASE09_'};

N=50;

for i = 1:length(speciesOrder)
    load(fullfile('..', 'outputs', [speciesOrder{i} '_exclude_not_found_and_G12'], 'ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1'))
    %[speciesGenomicData] = getGSMNMGenomeAnnotations(['MPI' speciesOrder{i} '-reactions_with_GPR.tsv']);
    %e = addGPRs(m.ensemble,speciesGenomicData.rxn_GPR_mapping, speciesFig{i});
    [speciesGenomicData] = getGSMNMGenomeAnnotations(['MPI' speciesOrder{i} '-reactions.tsv']);
    e = addGPRsGenBankIDs(m.ensemble,speciesGenomicData.rxn_GPR_mapping, speciesGenBank{i});
    ensembles(:,i) = e;
end

% ammonia;     urea;     L-Glutamic-Acid;L-Lysine; L-Serine
% 'cpd00013'; 'cpd00073'; 'cpd00023'; 'cpd00039'; 'cpd00054'
mediaCpdList = cellstr(['cpd00013'; 'cpd00023'; 'cpd00039'; 'cpd00054'; 'cpd00073']);
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
for k = 1:size(mediaCpdList,1);
  x = find(strcmp(seed_rxns_mat.mets, mediaCpdList(k,:)));
  if length(x) > 0
    mediaXSources = [mediaXSources, x];
  else
    fprintf(['AS-WARNING ' char(mediaCpdList(k,:)) ' (for media) not found in the rxn matrix\n']);
  end
end

mediaConditions = repmat(minimalMediaBase,[1,length(mediaXSources)]);
for i = 1:length(mediaXSources)
    mediaConditions(mediaXSources(i),i) = -5;
end


for i = 1:length(speciesOrder)
    curEnsemble = ensembles(:,i);
    
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
    
        fileName = ['../outputs/' speciesOrder{i} '_exclude_not_found_and_G12/geneEssentiality_' mediaCpdList{l} '.mat']; 
        save(fileName,'geneEssentialityByNet','allGenes', 'essentialGenes');
    end
end

