


speciesOrder = {'Root9', 'Root491', 'Root66D1'};
mediaCpdList = cellstr(['cpd00013'; 'cpd00023'; 'cpd00039'; 'cpd00054'; 'cpd00073']);
mediaNameList = cellstr(['ammo'; 'glut'; 'lysi'; 'seri'; 'urea']);

essgenes = zeros(3,5);
allgenes = zeros(3,5);
morethanquarter = zeros(3,5);
anygenes = zeros(3,5);
for i = 1:length(speciesOrder)
    for j = 1:length(mediaCpdList)
        load(fullfile('..', 'outputs', [speciesOrder{i} '_exclude_not_found_and_G12'], ['geneEssentiality_' mediaCpdList{j}]));
        %essGenes(i,j) = essentialGenes;
        %file = fopen(['../outputs/geneEssentiality/' speciesOrder{i} '_' mediaNameList{j} '.csv'],'w');
        %fprintf(file,'%s\n',essentialGenes{:});
        %fclose(file);
        essgenes(i,j) = size(essentialGenes, 1);
        allgenes(i,j) = size(allGenes, 1);
        N = sum(geneEssentialityByNet,2);
        anygenes(i,j) = size(allGenes(N > 0), 1);
        morethanquarter(i,j) = size(allGenes(N > 12), 1);
        T=array2table(N,'RowNames',allGenes);
        writetable(T,['../outputs/geneEssentiality/' speciesOrder{i} '_' mediaNameList{j} '.csv'],'WriteRowNames',true);
    end
end



