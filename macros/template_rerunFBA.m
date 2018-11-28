%************************************
%**  author: Antonella Succurro    **
%**  email:a.succurro[AT]gmail.com **
%**                                **
%**  created:       2018/06/04     **
%**  last modified: 2018/06/04     **
%************************************
%** Template to extract the biomass fluxes
%************************************
ensembleFname='XXXFNAME';
outpath='../outputs/XXXDNAME/';


% 'cpd00013'; 'cpd00073'; 'cpd00023'; 'cpd00039'; 'cpd00054'
myCpdList = cellstr([XXXCPDLIST]);

% Load the struct containing the ensemble
load(fullfile('..', 'outputs', 'XXXDNAME', ensembleFname))


load 2018_seed_rxns
seed_rxns_mat.X = -1*speye(length(seed_rxns_mat.mets));
seed_rxns_mat.Ex_names = strcat('Ex_',seed_rxns_mat.mets);

my_growth = computeFBAsol(m.ensemble, seed_rxns_mat, myCpdList, 'XXXCPDNAME', ensembleFname, outpath, 1, 1);
N=uint8(my_growth>0);
TG=array2table(N,'RowNames',myCpdList);
writetable(TG,[outpath ensembleFname, '_XXXCPDNAME_tab.csv'], 'WriteRowNames',true);

%N=biomass_fluxes;
%TG=array2table(N,'RowNames',allCpds);
%writetable(TG,[outpath ensembleFname, '_fba_growth.csv'], 'WriteRowNames',true);
