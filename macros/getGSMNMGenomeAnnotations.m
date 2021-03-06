function [GSMNMGenomicData] = getGSMNMGenomeAnnotations(genomeFile)
%--------------------------------------------------------------------------
% getGSMNMGenomeAnnotations - Reads in genome annotations, peg-to-rxn
% mappings, etc.
%
% Inputs:
%     universalRxnSet - Matlab structure containing an S matrix and a similar
%       matrix for exchange rxns (X matrix), a reversability indicator for
%       all rxns in S (rev), rxn IDs (rxns), rxn names (rxnNames), exchange
%       rxn names (Ex_names) metabolite IDs (mets), metabolite names (metNames),
%       and metabolite formulas (metFormulas)
%
% Outputs:
%     GSMNMGenomicData - a Matlab struct with the following fields:
%       rxn_GPR_mapping - a Matlab struct with the following fields:
%           rxns - a cell array of rxn IDs of the same size as "rxn_GPR_mapping.gprs"
%           gprs - a cell array of GPRs of the same size as "rxn_GPR_mapping.rxns"
%
% Written by Antonella Succurro, building on the function getPAGenomeAnnotation by Matt Biggs
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
% Import reactions from the draft reconstruction from KBase
%------------------------------------------------------------------------
% Relevant fields: id (reaction id); gpr (gene-protein-reaction)
rxns = tdfread(genomeFile);

%Keep only the first 8 chars (rxnXXXXX) NB biomass rxn is bio1
trimmed_rxnList = rxns.id(:,1:8);
trimmed_rxnList = strtrim(cellstr(trimmed_rxnList));

%gprfeat is the updated field with peg IDs from ModelSEED, as the KBase models had a different format
if isfield(rxns, 'gprfeat')
    trimmed_gprs = strtrim(cellstr(rxns.gprfeat));
else
    trimmed_gprs = strtrim(cellstr(rxns.gpr));
end

%Select only rxns with (non empty? check!) gpr field
%The only empty gpr field is the biomass, rxns w/o gpr have "Unknown"
keep = ~cellfun(@isempty,trimmed_gprs);
trimmed_rxnList = trimmed_rxnList(keep);
trimmed_gprs = trimmed_gprs(keep);

remove = find(strcmp(trimmed_gprs, 'Unknown'));
trimmed_rxnList(remove) = [];
trimmed_gprs(remove) = [];

rxn_GPR_mapping = struct;
rxn_GPR_mapping.rxns = trimmed_rxnList;
rxn_GPR_mapping.gprs = trimmed_gprs;

% Export
GSMNMGenomicData = struct;
GSMNMGenomicData.rxn_GPR_mapping = rxn_GPR_mapping;

end
