function [GSMNMData] = getGSMNMGrowthConditions(universalRxnSet, xSourcesTable, usePA14BM)
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
%     xSourcesTable - csv file for the growth matrix, containing the field SEEDID of the
%        Carbon/Nitrogen/Sulfur/Phosphate/X sources on which growth has been tested,
%        the Growth field (0/1 values) and the Well field relative to Biolog plates (optional)
%     usePA14BM - 0/1, if True the "standard" biomass composition will be used
% Inputs loaded from the environment:
%     biomass - function defining the biomass equation, obtained running the python script getBM.py
%     biomassFn_PA14 - standard biomass composition that was used in the original publication for all strains
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

% BIOMASS

switch nargin
  case 3
    usePA14BM = usePA14BM;
  otherwise
    usePA14BM = 1;
end

cpdsNotFound = [];

%Call the biomass function, automatically written by the script GetBM.py,
%to get the vectors of BM substrates and products
if usePA14BM
    load biomassFn_PA14_2018
    biomassFn = biomassFn_PA14;
else
    [bm_substrates, bm_products] = biomass();
    %The biomassFn is a vector with
    %+1 elements corresponding to products and
    %-1 elements corresponding to substrates
    %init as zeros:
    biomassFn = zeros(length(universalRxnSet.mets),1);

    % Match biomass' list of compound IDs to the seed_rxns_mat.mets indexes
    s = [];
    for k = 1:length(bm_substrates);
        x = find(strcmp(universalRxnSet.mets, bm_substrates(k)));
        if length(x) > 0;
            s = [s, x];
        else
            fprintf(['AS-WARNING ' char(bm_substrates(k)) ' biomass substrate not found!\n']);
            cpdsNotFound = [cpdsNotFound bm_substrates(k)];
        end
    end

    p = [];
    for k = 1:length(bm_products);
        x = find(strcmp(universalRxnSet.mets, bm_products(k)));
        if length(x) > 0
            p = [p, x];
        else
            fprintf(['AS-WARNING ' char(bm_products(k)) ' biomass product not found!\n']);
            cpdsNotFound = [cpdsNotFound bm_products(k)];
        end
    end

    biomassFn(s,1) = -1;
    biomassFn(p,1) = 1;
end


% MEDIA

%Call the baselinemedia function, automatically written by the script GetMedia.py,
%to get the vector metabolites composing the baseline media
blmedia = baselinemedia();

minimalMediaBase = zeros(length(universalRxnSet.mets),1);
b = [];
for k = 1:length(blmedia);
  b = [b, find(strcmp(universalRxnSet.mets, blmedia(k)))];
end

minimalMediaBase(b,1) = -1000;

%Open the xSourcesTable
xst = tdfread(xSourcesTable, ',');
if not (isfield(xst, 'SEEDID'))
  fprintf('\nERROR! The X source table does not contain the SEEDID field!\n\n')
end
if not (isfield(xst, 'Growth'))
  fprintf('\nERROR! The X source table does not contain the Growth field!\n\n')
end

% Get list of compound names sustaining growth
xGrowth = xst.SEEDID(xst.Growth > 0, :);
% Get indexes
growthXSources = [];
for k = 1:size(xGrowth,1);
  x = find(strcmp(universalRxnSet.mets, xGrowth(k,:)));
  if length(x) > 0
    growthXSources = [growthXSources, x];
  else
    fprintf(['AS-WARNING ' char(xGrowth(k,:)) ' (growth sustaining) not found in the rxn matrix\n']);
    cpdsNotFound = [cpdsNotFound cellstr(xGrowth(k,:))];
  end
end

% Get list of compound names to be included but not used for gapfilling
xNotForGapfill = xst.SEEDID(xst.Growth < 0, :);
% Get indexes
notForGapfillXSources = [];
for k = 1:size(xNotForGapfill,1);
  x = find(strcmp(universalRxnSet.mets, xNotForGapfill(k,:)));
  if length(x) > 0
    notForGapfillXSources = [notForGapfillXSources, x];
  else
    fprintf(['AS-WARNING ' char(xNotForGapfill(k,:)) ' (not for gapfill) not found in the rxn matrix\n']);
    cpdsNotFound = [cpdsNotFound cellstr(xNotForGapfill(k,:))];
  end
end

% Get list of compound names not sustaining growth
xNonGrowth = xst.SEEDID(xst.Growth == 0, :);
% Get indexes
nonGrowthXSources = [];
for k = 1:size(xNonGrowth,1);
  x = find(strcmp(universalRxnSet.mets, xNonGrowth(k,:)));
  if length(x) > 0
    nonGrowthXSources = [nonGrowthXSources, x];
  else
    fprintf(['AS-WARNING ' char(xNonGrowth(k,:)) ' (not growth sustaining) not found in the rxn matrix\n']);
    cpdsNotFound = [cpdsNotFound cellstr(xNonGrowth(k,:))];
  end
end

n = length(growthXSources);
growthConditions = repmat(minimalMediaBase,[1,n]);
for i = 1:n
    growthConditions(growthXSources(i),i) = -100;
end

nonGrowthConditions = repmat(minimalMediaBase,[1,length(nonGrowthXSources)]);
for i = 1:length(nonGrowthXSources)
    nonGrowthConditions(nonGrowthXSources(i),i) = -100;
end

notForGapfillConditions = repmat(minimalMediaBase,[1,length(notForGapfillXSources)]);
for i = 1:length(notForGapfillXSources)
    notForGapfillConditions(notForGapfillXSources(i),i) = -100;
end

GSMNMData = struct;
GSMNMData.biomassFn = biomassFn;
GSMNMData.minimalMediaBase = minimalMediaBase;
GSMNMData.growthXSources = growthXSources;
GSMNMData.growthConditions = growthConditions;
GSMNMData.nonGrowthXSources = nonGrowthXSources;
GSMNMData.nonGrowthConditions = nonGrowthConditions;
GSMNMData.notForGapfillXSources = notForGapfillXSources;
GSMNMData.notForGapfillConditions = notForGapfillConditions;

if length(cpdsNotFound) > 0
  printGrepLoop(cpdsNotFound);
end

end
