function [] = printGrepLoop(cellarray)

fprintf('for c in '); fprintf('%s ', cellarray{1:end}); fprintf('\n');
fprintf('do\ngrep $c simple-cpd-DB.tsv\ndone\n');
