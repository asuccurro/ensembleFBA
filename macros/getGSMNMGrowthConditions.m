function [GSMNMData] = getGSMNMGrowthConditions(universalRxnSet, xSourcesTable)
%--------------------------------------------------------------------------
% getGSMNMGrowthConditions - Generates a biomass function, and formats the
% growth conditions for use with the input rxn database.
%
% Inputs passed to the function:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames),
%       and metabolite formulas (metFormulas)
%     xSourcesTable - tsv file containing the SEEDID of the Carbon/Nitrogen/Sulfur/Phosphate/X
%        sources on which growth has been tested
% Inputs loaded from the environment:
%     biomass - function defining the biomass equation, obtained running the python script getBM.py
%     baselinemedia - function defining the baseline composition of the growth media, obtained
%     running the python script getMedia.py
% Outputs:
%     GSMNMData - a Matlab struct with the following fields:
%       biomassFn - the same format as a rxn (column) in universalRxnSet.S
%       growthXSources - the list of carbon/nitrogen/x sources which allow growth
%       growthConditions - a matrix of lower bounds for the exchange rxns in universalRxnSet.X
%       nonGrowthXSources - the list of carbon/nitrogen/x sources which do not allow growth
%       nonGrowthConditions - a matrix of lower bounds for the exchange rxns in universalRxnSet.X
%
% Written by Antonella Succurro, building on the function getPA14GrowthConditions by Matt Biggs
%--------------------------------------------------------------------------

% check that the Ex_names field exhists:

if not (isfield(universalRxnSet, 'Ex_names'))
  seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);
end

% BIOMASS

%Call the biomass function, automatically written by the script GetBM.py,
%to get the vectors of BM substrates and products
[bm_substrates, bm_products] = biomass();

%The biomassFn is a vector with
%+1 elements corresponding to products and
%-1 elements corresponding to substrates
%init as zeros:
biomassFn = zeros(length(universalRxnSet.Ex_names),1);

% Match biomass' list of compound IDs to the seed_rxns_mat.mets indexes
s = []
for k = 1:length(bm_substrates);
  s = [s, find(strcmp(universalRxnSet.mets, bm_substrates(k)))];
end

p = []
for k = 1:length(bm_products);
  p = [p, find(strcmp(universalRxnSet.mets, bm_products(k)))];
end

biomassFn(s,1) = -1;
biomassFn(p,1) = 1;

% MEDIA

%Call the baselinemedia function, automatically written by the script GetMedia.py,
%to get the vector metabolites composing the baseline media
blmedia = baselinemedia();

minimalMediaBase = zeros(length(universalRxnSet.Ex_names),1);
b = []
for k = 1:length(blmedia);
  b = [b, find(strcmp(universalRxnSet.mets, blmedia(k)))];
end

minimalMediaBase(b,1) = -1000;

%Open the xSourcesTable
xst = tdfread(xSourcesTable);
if not (isfield(xst, 'SEEDID'))
  fprintf('\nERROR! The X source table does not contain the SEEDID field!\n\n')
end
if not (isfield(xst, 'Growth'))
  fprintf('\nERROR! The X source table does not contain the Growth field!\n\n')
end

xGrowth = x.SEEDID
xNonGrowth = x.SEEDID
growthXSources = [];
for k = 1:length(xGrowth);
  growthXSources = [growthXSources, find(strcmp(universalRxnSet.mets, xGrowth(k,:)))];
end

n = length(growthXSources);
growthConditions = repmat(minimalMediaBase,[1,n]);
for i = 1:n
    growthConditions(growthXSources(i),i) = -10;
end

% Define non-growth conditions
% non-Growth X sources:
% 2,3-Butanediol        cpd01949    8832
% ACTN                  cpd00361    462
% CELB                  cpd00158    234
% fructose-6-phosphate  cpd00072    1859
% D-Galacturonate       cpd00280    8289
% Glucose-1-phosphate   cpd00089    6518
% glucose-6-phosphate   cpd00079    1866
% D-Glucarate           cpd00609    6445
% D-Mannose             cpd00138    1351
% Glucuronate           cpd00164    5261
% Sorbitol              cpd00588    7827
% Tartrate              cpd00666    4395
% Xylose                cpd00154    232
% Formate               cpd00047    7229
% Glycogen              cpd00155    231
% Glycolate             cpd00139    1350
% Glyoxalate            cpd00040    7224
% gly-asp-L             cpd11589    3438
% dAMP                  cpd00294    3272
% L-Methionine          cpd00060    6861
% Thyminose             cpd01242    3171
% Thymidine             cpd00184    4229
% TRHL                  cpd00794    7607
% L-Valine              cpd00156    230
% Maltose               cpd00179    602
% Amylotriose           cpd01262    3489
% Hydroxyphenylacetate  cpd03320    519
% 2-Oxobutyrate         cpd00094    1523
% L-Inositol            cpd00121    6001
% D-Mucic acid          cpd00652    8768
% L-Arabinose           cpd00224    7238
% L-Homoserine          cpd00227    7239
% SALC                  cpd00599    3162
% Citraconate           cpd01502    3003
% Acetoacetate          cpd00142    4919
% D-Malate              cpd00386    6024
% D-Ribose              cpd00105    5620
% D-Serine              cpd00550    2431
% Inosine               cpd00246    6219
% L-Threonine           cpd00161    5266

nonGrowthXSources = [8832,462,234,1859,8289,6518,1866,6445,1351,5261, ...
                          7827,4395,232,7229,231,1350,7224,3438,3272,6861, ...
                          3171,4229,7607,230,602,3489,519,1523,6001,8768, ...
                          7238,7239,3162,3003,4919,6024,5620,2431,6219,5266];

nonGrowthConditions = repmat(minimalMediaBase,[1,length(nonGrowthXSources)]);
for i = 1:length(nonGrowthXSources)
    nonGrowthConditions(nonGrowthXSources(i),i) = -10;
end

PA14Data = struct;
PA14Data.biomassFn = biomassFn;
PA14Data.minimalMediaBase = minimalMediaBase;
PA14Data.growthXSources = growthXSources;
PA14Data.growthConditions = growthConditions;
PA14Data.nonGrowthXSources = nonGrowthXSources;
PA14Data.nonGrowthConditions = nonGrowthConditions;

end
