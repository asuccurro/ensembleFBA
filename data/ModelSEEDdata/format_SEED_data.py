#------------------------------------------------------------------
# Read in SEED rxn database and format for easy use in MATLAB
# Written by Matt Biggs, 2016
# Modified by Antonella Succurro, 2018
#------------------------------------------------------------------

#------------------------------------------------------------------
# Classes and functions
#------------------------------------------------------------------
cmpdID2name = {}
cmpdID2formula = {}
rxnID2name = {}
usedCmpds = []
rxnIDlist = []

EXCLUDETRANSPORT = True

# original ModelSeed db files
#cmpdFileName='compounds.master_30Mar16.tsv'
#rxnFileName = 'reactions.master_30Mar16.tsv'

# updated ModelSeed db files
cmpdFileName='compounds.tsv'
rxnFileName='reactions.tsv'

#------------------------------------------------------------------
# Read compound file and make dictionary
#------------------------------------------------------------------
cmpdFile = open(cmpdFileName,'r')

for line in cmpdFile:
	if line.find('id	abbreviation	name	formula	mass	source') == -1:
		# 0 	1				2 		3		4		5		6			7		8		9			10				11			12		13			14	15	16					17				18
		# id	abbreviation	name	formula	mass	source	structure	charge	is_core	is_obsolete	linked_compound	is_cofactor	deltag	deltagerr	pka	pkb	abstract_compound	comprised_of	aliases
		lps = line.split('\t')
		cid = lps[0]
		# names
		name = lps[2]
		cmpdID2name[cid] = name
		# formula
		formula = lps[3]
		cmpdID2formula[cid] = formula
cmpdFile.close()

#------------------------------------------------------------------
# Make reaction network
#------------------------------------------------------------------
rxnFile = open(rxnFileName,'r')

# First, keep only the compounds from valid reactions
usedCmpds = []
for line in rxnFile:
	if line.find('id	abbreviation	name	code	stoichiometry	is_transport	equation') == -1:
		# 0		1				2		3		4				5				6			7			8				9			10					11			12		13			14		15			16				17		18			19
		# id	abbreviation	name	code	stoichiometry	is_transport	equation	definition	reversibility	direction	abstract_reaction	pathways	aliases	ec_numbers	deltag	deltagerr	compound_ids	status	is_obsolete	linked_reaction
		lps = line.split('\t')
		is_transport = lps[5].strip()
		status = lps[17].strip()
		is_obsolete = lps[18].strip()
		if status.find('OK') > -1 and \
		   status.find('CI') == -1 and \
		   status.find('MI') == -1 and \
		   is_obsolete == '0' and (is_transport == '0' and EXCLUDETRANSPORT):
			rxnID = lps[0]
			rxnName = lps[2]
			rxnID2name[rxnID] = rxnName
			rxnIDlist.append(rxnID)
			stoichiometry = lps[4]
			# Example: -1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"
			stoich_parts = stoichiometry.split(';')
			for sp in stoich_parts:
				spps = sp.split(':')
				usedCmpds.append(spps[1])
rxnFile.close()

usedCmpds = list(set(usedCmpds))
rxnIDset = set(rxnIDlist)

if EXCLUDETRANSPORT:
    print 'Excluding transport or obsolete reactions, in the complete database there are:'
else:
    print 'NOT Excluding transport, excluding obsolete reactions, in the complete database there are:'
print str(len(rxnIDlist)) + ' reactions'
print str(len(usedCmpds)) + ' compounds'

# Next, write reactions to file in a sparse format
rxnFile = open(rxnFileName,'r')
pathMatOutFile = open('complete_SEED_matrix.tsv','w')

for line in rxnFile:
	if line.find('id	abbreviation	name	code	stoichiometry	is_transport	equation') == -1:
		# 0		1				2		3		4				5				6			7			8				9			10					11			12		13			14		15			16				17		18			19
		# id	abbreviation	name	code	stoichiometry	is_transport	equation	definition	reversibility	direction	abstract_reaction	pathways	aliases	ec_numbers	deltag	deltagerr	compound_ids	status	is_obsolete	linked_reaction
		lps = line.split('\t')
		rxnID = lps[0]
		status = lps[17].strip()
		is_obsolete = lps[18].strip()
		if rxnID in rxnIDset:
			stoichiometry = lps[4]

			rev = '0'
			if lps[8].find('='):
				rev = '1'

			# Write reaction in sparse format
			stoich_parts = stoichiometry.split(';')
			sparseStr = rxnID + '\t' + rev + '\t'
			for sp in stoich_parts:
				spps = sp.split(':')
				# each pair in parentheses contains the index of the coefficient, and the value at that index
				sparseStr += '(' + str(usedCmpds.index(spps[1])+1) + ',' + spps[0] + ')|'
			sparseStr = sparseStr.rstrip('|')
			sparseStr += '\t' + rxnID2name[rxnID]
			pathMatOutFile.write(sparseStr + '\n')

rxnFile.close()
pathMatOutFile.close()

# Other output files
cmpdNameOutfile = open('compound_info_SEED.tsv','w')
for cmpID in usedCmpds:
	cmpdNameOutfile.write(cmpID + '\t' + cmpdID2formula[cmpID] + '\t' + cmpdID2name[cmpID] + '\n')
cmpdNameOutfile.close()
