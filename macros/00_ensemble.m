%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/05/24     **
%**  last modified: 2018/06/04     **
%************************************
% Set parameters
params = struct;
params.ISTEST = XXXTEST;
params.stochast = XXXSTOC;
params.numModels2gen = XXXSIZE;
params.verbose = 0;

params.fileGrowthConditions='growthMatrix_Root491_XXXCOND.csv';
params.fileOutPath='../outputs/root491/XXXCOND/'

params.fileAnnotations='rhizobiumRoot491-reactions.tsv';
params.sequential = 1;
params.numModels2gen = 1;
params.fractionUrxns2set = 0.8;
params.rndSequence = 1;

%run ensemble
[ens] = runEnsemble(params);
