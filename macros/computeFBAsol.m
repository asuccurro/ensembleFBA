function [fba_growth] = computeFBAsol(ensemble, rxns_mat, fbaCpdList, fbaCpdName, ensembleFname, outpath, savesol, isminimalmedia)
% Written by Antonella Succurro, building on the function getPA14GrowthConditions by Matt Biggs
%--------------------------------------------------------------------------
fbaXSources = [];
for k = 1:size(fbaCpdList,1);
  x = find(strcmp(rxns_mat.mets, fbaCpdList(k,:)));
  if length(x) > 0
    fbaXSources = [fbaXSources, x];
  else
    fprintf(['AS-WARNING ' char(fbaCpdList(k,:)) ' (for fba) not found in the rxn matrix\n']);
  end
end

if isminimalmedia
    fprintf('Using minimal media formulation\n');
    blmedia = minimalmedia();
else
    blmedia = baselinemedia();
end

minimalMediaBase = zeros(length(rxns_mat.mets),1);
b = [];
for k = 1:length(blmedia);
  b = [b, find(strcmp(rxns_mat.mets, blmedia(k)))];
end

minimalMediaBase(b,1) = -100;

% Manual fixes
% CO2
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00011')), 1) = 0;
% Limiting nutrients
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00009')), 1) = -5;
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00027')), 1) = -5;
minimalMediaBase(find(strcmp(rxns_mat.mets,  'cpd00048')), 1) = -5;

fbaConditions = repmat(minimalMediaBase,[1,length(fbaXSources)]);
for i = 1:length(fbaXSources)
    fbaConditions(fbaXSources(i),i) = -5;
end

[fba_growth] = ensembleFBA(ensemble,rxns_mat.Ex_names,fbaConditions,0);
N=fba_growth;
TG=array2table(N,'RowNames',fbaCpdList);
writetable(TG,[outpath ensembleFname, '_fba_', fbaCpdName, '_growth.csv'], 'WriteRowNames',true);

if savesol
    
    e = ensemble;
    r = rxns_mat.Ex_names;
    c = fbaConditions;

    % for every media condition and for every network store the solution fluxes
    % store in place 1 the reactions, in place 2 the EX_rxns
    solutions = cell(size(c,2),2);
    for i = 1:size(c,2)
        s_rx = zeros(size(rxns_mat.rxns,1),length(e));
        s_ex = zeros(size(rxns_mat.Ex_names,1),length(e));
        for j = 1:length(e)
            % match model reactions to matrix; distinguish exchange from rxns as are named differently
            model = e{j};
            rxnindex = [];
            exrindex = [];
            rr = [];
            ee = [];
            for z = 1:length(model.rxns)
                ef = find(strcmp(rxns_mat.Ex_names, model.rxns(z)));
                rf = find(strcmp(rxns_mat.rxns, model.rxns(z)));
                if length(rf) > 0
                    rxnindex = [rxnindex rf];
                    rr = [rr z];
                elseif length(ef)>0
                    exrindex = [exrindex ef];
                    ee = [ee z];
                end
            end
            [growth,x] = fba_flex(model,r,c(:,i),1);
            s_rx(rxnindex,j) = x(rr);
            s_ex(exrindex,j) = x(ee);
        end
        solutions{i,1} = s_rx;
        solutions{i,2} = s_ex;
        % Actually we do not care about Ex reactions
        N = s_rx;
        T=array2table(N,'RowNames',rxns_mat.rxns);
        writetable(T,[outpath ensembleFname '_fba_sol_' char(fbaCpdList(i)) '.csv'], 'WriteRowNames',true);
        dlmwrite([outpath ensembleFname, '_exc_rxns_' char(fbaCpdList(i)) '.csv'], s_ex, ',');
    end
end
