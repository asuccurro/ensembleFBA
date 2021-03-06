function [newModel,exchangesNotAdded] = addExchangeRxns(model,metList)
%-------------------------------------------------------------------------- 
% addExchangeRxns - Add exchange reactions to a model using a list of SEED 
% metabolite IDs.
%
% Inputs:
%     model - COBRA-format metabolic network model (Matlab struct)
%     metList - a cell array of SEED metabolite identifiers
%
% Outputs:
%     newModel - COBRA-format metabolic network model (Matlab struct)
%     exchangesNotAdded - list of metabolites that cannot be added
%
% Written by Matt Biggs, mb3ad@virginia.edu, 2016
%-------------------------------------------------------------------------- 

metsInModel = metList(ismember(metList,model.mets));
exchangesNotAdded = metList(~ismember(metList,model.mets));

% Add exchange reactions
newModel = model;
for i = 1:length(metsInModel)
    % Check if exchanges already exist
    exRxnIndex = find(ismember(model.rxns,['Ex_' metsInModel{i}]));
    
    % If exchange does not already exist, add it
    if isempty(exRxnIndex)
        curMetIndex = find(ismember(model.mets,metsInModel{i}));
        tmpExRxn = zeros(size(model.mets));
        tmpExRxn(curMetIndex) = -1;
        newModel.S(:,end+1) = tmpExRxn;
        newModel.rxns{end+1} = ['Ex_' metsInModel{i}];
        newModel.rxnNames{end+1} = ['Ex_' metsInModel{i}];
        newModel.rev(end+1) = 1;
        newModel.c(end+1) = 0;
        newModel.ub(end+1) = 0;
        newModel.lb(end+1) = 0;
        newModel.grRules{end+1} = '';
    end
end


end