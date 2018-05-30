function [biomassFn] = myBiomassBuilder(universalRxnSet, bm_substrates, bm_products)

biomassFn = zeros(length(universalRxnSet.mets),1);

% Match biomass' list of compound IDs to the universalRxnSet.mets indexes
s = [];
for k = 1:length(bm_substrates);
  x = find(strcmp(universalRxnSet.mets, bm_substrates(k)));
  if length(x) > 0;
    s = [s, x];
  else
    fprintf(['AS-WARNING ' char(bm_substrates(k)) ' biomass substrate not found!\n']);
  end
end

p = [];
for k = 1:length(bm_products);
  x = find(strcmp(universalRxnSet.mets, bm_products(k)));
  if length(x) > 0
    p = [p, x];
  else
    fprintf(['AS-WARNING ' char(bm_products(k)) ' biomass product not found!\n']);
  end
end

biomassFn(s,1) = -1;
biomassFn(p,1) = 1;
