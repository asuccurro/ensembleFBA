%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
%** Template to extract the biomass fluxes
%************************************
ensembleFname='ensemble_50_size_26_gcs_11_ngcs_stochasticWeights_1';
outpath='../outputs/test/';


% 'cpd00013'; 'cpd00073'; 'cpd00023'; 'cpd00039'; 'cpd00054'
myCpdList = cellstr(['cpd00013']);

% Load the struct containing the ensemble
load(fullfile('..', 'outputs', 'test', ensembleFname))


load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

my_growth = computeFBAsol(m.ensemble, seed_rxns_mat, myCpdList, 'test', ensembleFname, outpath, 1, 1);
N=uint8(my_growth>0);
TG=array2table(N,'RowNames',myCpdList);
writetable(TG,[outpath ensembleFname, '_minimal_tab.csv'], 'WriteRowNames',true);


my_growth = computeFBAsol(m.ensemble, seed_rxns_mat, myCpdList, 'test', ensembleFname, outpath, 1, 0);
N=uint8(my_growth>0);
TG=array2table(N,'RowNames',myCpdList);
writetable(TG,[outpath ensembleFname, '_biolog_tab.csv'], 'WriteRowNames',true);

%N=biomass_fluxes;
%TG=array2table(N,'RowNames',allCpds);
%writetable(TG,[outpath ensembleFname, '_fba_growth.csv'], 'WriteRowNames',true);
