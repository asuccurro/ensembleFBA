# Genome Scale Model Analysis of MPI Root9, Root491, and Root66D1 with EnsembleFBA

This is a forked copy of Biggs MB and Papin JA's EnsembleFBA code, adapted for the
analysis performed in [Jacoby, Succurro, Kopriva 2018](MANUSCRIPT). For the original EnsembleFBA repository [click here](https://github.com/mbi2gs/ensembleFBA).
For questions contact asuccurro (at) protonmail (dotcom).


## Setup

EnsembleFBA scripts are either written in MATLAB, Python or bash. You will therefore need MATLAB license etc.

Once cloned, export the local path, e.g.
```bash
export ensembleFBA=/home/succurro/repositories/publicgithub/ensembleFBA/
```


### 1. Make sure gurobi is installed (optimization problems) and linked to MATLAB

1. [Download](http://www.gurobi.com/downloads/gurobi-optimizer)
2. [Install](http://www.gurobi.com/documentation/8.0/quickstart_linux/software_installation_guid.html)
3. [Get a license](http://www.gurobi.com/downloads/licenses/license-center)

```bash
cd $GUROBI_HOME/matlab
matlab -nodesktop -nosplash -r "gurobi_setup"
```

Remember to add the path $GUROBI_HOME in matlab (see `$ensembleFBA/runmatlab/startup.m` file)

### 2. Setup a virtualenv for the python scripts

```bash
cd $ensembleFBA
virtualenv -p /usr/bin/python3 venvpy
source venvpy/bin/activate
pip install -r req_venvpy.txt
```

### 3. Test runs


```bash
matlab -nodesktop -nosplash -nodisplay  -r "testMPIRoot('Root9'); quit"
matlab -nodesktop -nosplash -nodisplay  -r "testMPIRoot('Root491'); quit"
matlab -nodesktop -nosplash -nodisplay  -r "testMPIRoot('Root66D1'); quit"
```


### 4. Other technical details

The above setup are enough to re-run the simulations of [Jacoby, Succurro, Kopriva 2018](MANUSCRIPT).
The additional instructions [at the end of the wiki](#instructions-for-novel-runs) are meant to enhance reproducibility and to run new/updated simulation.
Content:
1. [Updating the ModelSEED db](#modelseed-db)
2. [Defintion of a biomass function](#biomass-function)
3. [Reconstruction of draft networks with KBase](#draft-networks)
4. [Definition of growth matrices from biolog plate experiments](#growth-matrices)
5. [Test reproducibility with fixed random seed](#test-reproducibility)

## EnsembleFBA workflow

Excluding the [compounds not found in the matrix and G12](#note-on-n-sources), the numbers of growth and non growth media are:

| |  Root9 | Root491 | Root66D1 |
| --- | --- | --- | --- |
| 1 | 54 | 58 | 52 |
| 0 | 26 | 22 | 28 |


### Ensemble generation

I made a conservative choice for the number of training conditions (gapfilling media):
half of the lowest number of growth-sustaining or growth-sustaining compounds, i.e. 26 and 11 respectively.

The ensembles are generated in two versions: excluding not found compounds and G12 (`growthMatrix_Root9_exclude_not_found_and_G12.csv`),
and excluding also the 5 N sources with proteomics data (`growthMatrix_Root9_exclude_5N_and_nf_and_G12.csv`)

Running:
```bash
source runEnsemble_manuscript.sh
```
will start in total 6 reconstruction jobs in the background. On a 16 core desktop, they take about 2 days to complete.
The jobs can be monitored via the log files (`/tmp/ensemble_Root*`) and error files (`/tmp/err_ensemble_Root*`).


### Analysis on ensembles

Once the ensembles are successfully generated, ensembleFBA will perform FBA on all the networks on all the conditions.
The script:

```bash
source analyseEnsemble_manuscript.sh
```

will produce various files, all named starting with "ensemble_X_size_Y_gcs_Z_ngcs_stochasticWeights_1". The files are described below.

#### _conditions.csv



#### _gc_tab.csv; _ngc_tab.csv; _nfg_tab.csv


#### _proteomics_growth.csv


#### _fba_sol_cpdXXXXX.csv and _exc_rxns_cpdXXXXX.csv











## Instructions for novel runs

### ModelSEED db

The database for compounds and reactions were downloaded from github.com/ModelSEED/ModelSEEDDatabase/ on May 24th, 2018 from
the master branch (corresponding in the current rebased master to commit 25c8d32e50cefe16f3bee5725577c6035ca0a5b8), and the reaction matrix updated as follows:

```bash
cd $ensembleFBA/data/ModelSEEDdata
# Obtain the latest files from ModelSEED
wget https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.tsv
wget https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/reactions.tsv
# Run format_SEED_data.py to obtain complete_SEED_matrix.tsv and compound_info_SEED.tsv
python format_SEED_data.py
cd $ensembleFBA/runmatlab
# Run package_SEED_data.m to obtain 2018_seed_rxns.mat (moved then to $ensembleFBA/data/ModelSEEDdata)
source generateSEEDMat.sh
# Done! :)
```

### Biomass function

Following what done in the original EnsembleFBA analysis, a standard (minimal) biomass function definition is used. The
biomass matrix (file `$ensembleFBA/data/MPIRoots/biomassFn_PA14_2018.mat`) has to
be updated with the new database, and is obtained running:

```bash
cd $ensembleFBA/runmatlab
source updateBM.sh
```

### Draft networks

Draft, non-gapfilled Genome Scale Metabolic Network Models were obtained through [this KBase narrative](https://narrative.kbase.us/narrative/ws.37070.obj.1).
The tsv files were downloaded and re-organized into the folders

```bash
data/MPIRoot491/MPIRoot491-compounds.tsv
data/MPIRoot491/MPIRoot491-reactions.tsv
data/MPIRoot9/MPIRoot9-compounds.tsv
data/MPIRoot9/MPIRoot9-reactions.tsv
data/MPIRoot66D1/MPIRoot66D1-compounds.tsv
data/MPIRoot66D1/MPIRoot66D1-reactions.tsv
```

These files needed re-formatting to solve a problem with delimiters in MATLAB:

```bash
cd $ensembleFBA/macros
./reformatReactions.py -i Root9
./reformatReactions.py -i Root491
./reformatReactions.py -i Root66D1
```

### Growth matrices

Define growth/no growth media from input experimental data (see `data/MPIRoots/biolog_summary.tsv`). In this case, the data are from a Biolog plate with
95 different N sources. The media composition is taken from the file [`media_list_with_meta.tsv`](https://github.com/ModelSEED/ModelSEEDDatabase/blob/master/Media/media_list_with_meta.txt)
selecting the "Nitrogen" media. The script `produceGrowthMat.py` runs different growth matrix constructions.

```bash
cd $ensembleFBA/macros
# Write the matlab file defining the baseline media and the table of N sources
./getMedia.py -m '../data/ModelSEEDdata/media_list_with_meta.tsv' -c 'Nitrogen' -b '../data/MPIRoots/biolog_summary.tsv' -o '../data/MPIRoots/singleNMedia/' -d '../data/ModelSEEDdata/compounds.tsv' -F
# Write the growth matrix for the selected organism on the selected N sources
python produceGrowthMat.py
# Done :)
```

Output files:
* in data/MPIRoots/: individual minimal media with different N source files; `ncompounds.tsv` table with the identified N compounds.
* in the respective organism folder: `growthMatrix_RootX.csv` table with growth/no growth per N source


#### Note on N sources

Out of the 95 N sources, 84 are automatically matched to the DB.
F7 (Guanosine) was removed for experimental issues, and 
10 other are not found by the `GetMedia.py` macro in the file `data/ModelSEEDdata/compounds.tsv`.
After checking manually, they could be matched.
However, out of these 10 compounds, 9 are then not found in the reaction matrix, as a consequence of the exclusion of transport reactions.
This might suggest annotation issues for these compounds. The only compound included, G12 (Valeric Acid),
is then found to cause infeasibility when building the ensemble. 4 more compounds (C4, C9, G7 and E10) are also not found in the reaction matrix.
The final choice is to exclude these 14 compounds:

| Well | Name | Class | ID |
| --- | --- | --- | --- |
| C4 | D-Asparagine | D-amino acids | cpd01308 |
| C9 | D-Valine | D-amino acids | cpd03840 |
| D2 | N-Phthaloyl-L-Glutamic-Acid | Amino acid derivatives | cpd24431 |
| D6 | N-Amylamine | Amines | cpd24432 |
| D7 | N-Butylamine | Amines | cpd19969 |
| D10 | Ethylenediamine | Amines | cpd24433 |
| E6 | Glucuronamide | Modified sugars | cpd26269 |
| E7 | D-L-Lactamide | Modified sugars | cpd23860 |
| E10 | D-Mannosamine | Modified sugars | cpd02241 |
| G4 | Alloxan | Nitrogen bases | cpd24434 |
| G6 | Parabanic-Acid | Other organic compounds | cpd24435 |
| G7 | D-L-a-Amino-N-Butyric-Acid | Amino acid derivatives | cpd01573 |
| G10 | D-L-a-Amino-Caprylic-Acid | Amino acid derivatives | cpd24436 |
| G12 | a-Amino-N-Valeric-Acid | Amino acid derivatives | cpd01258 |


### Test reproducibility

Run the `reproduce_runEnsemble_allRoots.sh` script, which sets the random seed and has verbosity turned on.
Check that the networks are reconstructed on the same conditions in the same order:

```bash
for f in /tmp/ensemble_Root66D1_*10_21*; do awk '/p1/{n=NR} n && NR==n+4 || NR==n+8' $f > /tmp/gcorder_${f:5:-3}out; done
```

will print the conditions in separate files. Those files should be identical for the conditions `_exclude_G12` and `_exclude_not_found_and_G12`, while the condition not excluding any compound should at some point differ (the condition "55" should at some point appear).


