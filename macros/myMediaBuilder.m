function [minimalMediaBase, growthConditions, nonGrowthConditions, growthXSources, nonGrowthXSources] = myMediaBuilder(universalRxnSet, blmedia, xGrowth, xNonGrowth)

minimalMediaBase = zeros(length(universalRxnSet.mets),1);
b = [];
for k = 1:length(blmedia);
  b = [b, find(strcmp(universalRxnSet.mets, blmedia(k)))];
end

minimalMediaBase(b,1) = -1000;


growthXSources = [];
for k = 1:length(xGrowth);
  x = find(strcmp(universalRxnSet.mets, xGrowth(k,:)));
  if length(x) > 0
    growthXSources = [growthXSources, x];
  else
    fprintf(['AS-WARNING ' char(xGrowth(k,:)) ' (growth sustaining) not found in the rxn matrix\n']);
  end
end

nonGrowthXSources = [];
for k = 1:length(xNonGrowth);
  x = find(strcmp(universalRxnSet.mets, xNonGrowth(k,:)));
  if length(x) > 0
    nonGrowthXSources = [nonGrowthXSources, x];
  else
    fprintf(['AS-WARNING ' char(xNonGrowth(k,:)) ' (not growth sustaining) not found in the rxn matrix\n']);
  end
end

n = length(growthXSources);
growthConditions = repmat(minimalMediaBase,[1,n]);
for i = 1:n
    growthConditions(growthXSources(i),i) = -10;
end

nonGrowthConditions = repmat(minimalMediaBase,[1,length(nonGrowthXSources)]);
for i = 1:length(nonGrowthXSources)
    nonGrowthConditions(nonGrowthXSources(i),i) = -10;
end
