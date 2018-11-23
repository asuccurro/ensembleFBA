import getGrowthMat


biolog    = '../data/MPIRoots/biolog_summary.tsv'
compounds = '../data/MPIRoots/singleNMedia/ncompounds.tsv'
excl      = ''
selc      = ''
for orgid in ['Root9', 'Root491', 'Root66D1']:
    outpath = '../data/MPI'+orgid+'/'
    for classid in ['Modified sugars', 'Oligopeptides', 'Other organic compounds', 'Amines', 'D-amino acids', 'Inorganic compounds', 'Amino acid derivatives', 'L-amino acids', 'Nucleosides and nucleotides', 'Nitrogen bases']:
        getGrowthMat.getGrowthMatrix(biolog, compounds, outpath, orgid, classid, excl, selc)
    getGrowthMat.getGrowthMatrix(biolog, compounds, outpath, orgid, 'NA', 'A2 A5 A12 B6 B10', selc)
    getGrowthMat.getGrowthMatrix(biolog, compounds, outpath, orgid, 'NA', excl, selc)
    getGrowthMat.getGrowthMatrix(biolog, compounds, outpath, orgid, 'NA', 'G12', selc)
    getGrowthMat.getGrowthMatrix(biolog, compounds, outpath, orgid, 'NA', 'C4 C9 D2 D6 D7 D10 E6 E7 E10 G4 G6 G7 G10 G12', selc)
    getGrowthMat.getGrowthMatrix(biolog, compounds, outpath, orgid, 'NA', 'A2 A5 A12 B6 B10 C4 C9 D2 D6 D7 D10 E6 E7 E10 G4 G6 G7 G10 G12', selc)

