
load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);



% Get the Streptococcus data formatted to work with the SEED database
[StrepData] = getStrepGrowthConditions(seed_rxns_mat);

% Get the Streptococcus gene/drug mappings
[StrepGeneDrugPairs] = getStrepGeneDrugPairs();


% load the ensembles of the 3 strains
ensembles = cell(50,3);
speciesOrder = {'Root9', 'Root491', 'Root66D1'};

for i = 1:length(speciesOrder)
    load(fullfile('..', 'outputs', [speciesOrder{i} '_exclude_not_found_and_G12'], 'ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1'))
    [speciesGenomicData] = getGSMNMGenomeAnnotations(['MPI' speciesOrder{i} '-reactions.tsv']);
    e = addGPRs(m.ensemble,speciesGenomicData.rxn_GPR_mapping);
    ensembles(:,i) = e;
end

% ammonia; ; ; L-Lysine; 
% 'cpd00013'; 'cpd00073'; 'cpd00023'; 'cpd00039'; 'cpd00054'
mediaCpdList = cellstr(['cpd00013'; 'cpd00039']);
blmedia = minimalmedia();
minimalMediaBase = zeros(length(rxns_mat.mets),1);

x = [];
for k = 1:length(blmedia);
  x = [x, find(strcmp(rxns_mat.mets, blmedia(k)))];
end

minimalMediaBase(x,1) = -100;
% Limiting nutrients
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00009')), 1) = -5;
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00027')), 1) = -5;
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00048')), 1) = -5;

mediaXSources = [];
for k = 1:size(mediaCpdList,1);
  x = find(strcmp(rxns_mat.mets, mediaCpdList(k,:)));
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



for i = 1:length(rxnMappingsList)
    curEnsemble = ensembles(:,i);
    
    allGenes = cell(0,1);
    for j = 1:N
        allGenes = [allGenes; curEnsemble{j}.genes];
    end
    allGenes = unique(allGenes);
    
    geneEssentialityByNet = zeros(length(allGenes),N);
    for j = 1:N
        curMod = curEnsemble{j};
        for k = 1:length(allGenes)
            curGene = allGenes{k};
            delMod = simulateGeneDeletion(curMod,curGene);
            delGrowth = fba_flex(delMod,seed_rxns_mat.Ex_names,richMedia,0);
            geneEssentialityByNet(k,j) = delGrowth < 1e-10;
        end
    end
    
    % Decide on essential genes
    essentialGeneIndicators = sum(geneEssentialityByNet,2) > N/2;
    essentialGenes = allGenes(essentialGeneIndicators > 0);
    drugsHitEssentialGenes = genesWdrugMatches.drugIDs(ismember(genesWdrugMatches.genes,essentialGenes));
    
    fileName = ['CE13_geneEssentiality_' speciesOrder{i} '.mat']; fileName = strrep(fileName,'S. ','');
    save(fileName,'geneEssentialityByNet','allGenes','drugsHitEssentialGenes');
end

% Evaluate the drugs that uniquely interact with each species (and the
% genes they target)
load CE13_geneEssentiality_mitis.mat
mitis_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       mitis_drugs = [mitis_drugs; curDrug];
   end
end
mitis_drugs = unique(mitis_drugs);

load CE13_geneEssentiality_gallolyticus.mat
gallolyticus_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       gallolyticus_drugs = [gallolyticus_drugs; curDrug];
   end
end
gallolyticus_drugs = unique(gallolyticus_drugs);

load CE13_geneEssentiality_oralis.mat
oralis_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       oralis_drugs = [oralis_drugs; curDrug];
   end
end
oralis_drugs = unique(oralis_drugs);

load CE13_geneEssentiality_equinus.mat
equinus_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       equinus_drugs = [equinus_drugs; curDrug];
   end
end
equinus_drugs = unique(equinus_drugs);

load CE13_geneEssentiality_pneumoniae.mat
pneumoniae_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       pneumoniae_drugs = [pneumoniae_drugs; curDrug];
   end
end
pneumoniae_drugs = unique(pneumoniae_drugs);

load CE13_geneEssentiality_vestibularis.mat
vestibularis_drugs = {};
for i = 1:length(drugsHitEssentialGenes)
   curDrugList =  drugsHitEssentialGenes{i};
   splitDrugList = strsplit(curDrugList,',');
   for j = 1:length(splitDrugList)
       curDrug = splitDrugList{j};
       vestibularis_drugs = [vestibularis_drugs; curDrug];
   end
end
vestibularis_drugs = unique(vestibularis_drugs);

size(mitis_drugs)
size(gallolyticus_drugs)
size(oralis_drugs)
size(equinus_drugs)
size(pneumoniae_drugs)
size(vestibularis_drugs)

unique2mitis = mitis_drugs(~ismember(mitis_drugs,[gallolyticus_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2gallolyticus = gallolyticus_drugs(~ismember(gallolyticus_drugs,[mitis_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2oralis = oralis_drugs(~ismember(oralis_drugs,[mitis_drugs;gallolyticus_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2equinus = equinus_drugs(~ismember(equinus_drugs,[mitis_drugs;gallolyticus_drugs;oralis_drugs;pneumoniae_drugs;vestibularis_drugs]))
unique2pneumoniae = pneumoniae_drugs(~ismember(pneumoniae_drugs,[mitis_drugs;gallolyticus_drugs;oralis_drugs;equinus_drugs;vestibularis_drugs]))
unique2vestibularis = vestibularis_drugs(~ismember(vestibularis_drugs,[mitis_drugs;gallolyticus_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs]))

allDrugs = [mitis_drugs;gallolyticus_drugs;oralis_drugs;equinus_drugs;pneumoniae_drugs;vestibularis_drugs];
uDrugs = unique(allDrugs);
counts = zeros(size(uDrugs));
for i = 1:length(uDrugs)
    counts(i) = sum(ismember(allDrugs,uDrugs{i}));
end
common2all = uDrugs(counts == 6)

%---------------------------------------------------------------------
% Evaluate the essential reactions in each species and the subsystem
% enrichment in each
%---------------------------------------------------------------------
% S. mitis
inputFile = 'EssRxns_S.mitis_';
[mitis_subsysPvals,uSubsys] = calcStrepSubsystemEnrichment(inputFile,N);

subsysEnrichment = zeros(length(uSubsys),6);
subsysEnrichment(:,1) = mitis_subsysPvals;

fid = fopen('CE13_uniqueSubsystems.txt','w');
for i = 1:length(uSubsys)
   fprintf(fid,[uSubsys{i} '\n']); 
end
fclose(fid);

% S. gallolyticus
inputFile = 'EssRxns_S.gallolyticus_';
[gallolyticus_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,2) = gallolyticus_subsysPvals;

% S. oralis
inputFile = 'EssRxns_S.oralis_';
[oralis_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,3) = oralis_subsysPvals;

% S. equinus
inputFile = 'EssRxns_S.equinus_';
[equinus_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,4) = equinus_subsysPvals;

% S. pneumoniae
inputFile = 'EssRxns_S.pneumoniae_';
[pneumoniae_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,5) = pneumoniae_subsysPvals;

% S. vestibularis
inputFile = 'EssRxns_S.vestibularis_';
[vestibularis_subsysPvals,~] = calcStrepSubsystemEnrichment(inputFile,N);
subsysEnrichment(:,6) = vestibularis_subsysPvals;


