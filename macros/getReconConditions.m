%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
function [conditionList] = getReconConditions(universalRxnSet, network)
%--------------------------------------------------------------------------
%
% Inputs passed to the function:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames),
%       and metabolite formulas (metFormulas)
%     network - struct containing one element of the ensemble
% Outputs:
%     conditionList - struct of the conditions used for the reconstruction
%--------------------------------------------------------------------------

gci = ~(network.growthConditions==0);
gc = [];
for j=1:size(gci,2)
  gc = [gc setdiff(universalRxnSet.mets(gci(:,j)), universalRxnSet.mets(gci(:,1+mod(j,size(gci,2)))))];
end
ngci = ~(network.nonGrowthConditions==0);
ngc = [];
for j=1:size(ngci,2)
  ngc = [ngc setdiff(universalRxnSet.mets(ngci(:,j)), universalRxnSet.mets(ngci(:,1+mod(j,size(ngci,2)))))];
end

conditionList = struct;
conditionList.growth = gc;
conditionList.nonGrowth = ngc;
end
