load 2018_seed_rxns
load biomassFn_PA14_2018

S = seed_rxns_mat.mets(biomassFn_PA14 < 0);
P = seed_rxns_mat.mets(biomassFn_PA14 > 0);

dlmwrite('../data/MPIRoots/biomassFn_PA14_substrates.csv', S, '');
dlmwrite('../data/MPIRoots/biomassFn_PA14_products.csv', P, '');
