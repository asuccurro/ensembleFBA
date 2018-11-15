%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/05/24     **
%**  last modified: 2018/06/04     **
%************************************
% Set parameters
params = struct;

%XXXORG Root9
%XXXCOND _exclude_A2-A5-A12-B6-B10

params.ISTEST = XXXTEST;
params.stochast = XXXSTOC;
params.numModels2gen = XXXSIZE;
params.numGC = XXXGC;
params.numNGC = XXXNGC;

params.fileGrowthConditions='growthMatrix_XXXORGXXXCOND.csv';
params.fileOutPath='../outputs/XXXORGXXXCOND/';

params.fileAnnotations='MPIXXXORG-reactions.tsv';

params.verbose = 1;
params.rndSeed = 15112018;

params.sequential = 1;
params.fractionUrxns2set = 0.8;
params.rndSequence = 1;


%run ensemble
[ens] = runEnsemble(params);
