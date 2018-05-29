#!/usr/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/25     **
#**  last modified: 2018/05/25     **
#************************************

import argparse
import subprocess
import collections

def main():
    args = options()
    verbose = args.verbose
    # sbp = subprocess.Popen(("grep", "-w", "%s"%args.bmname, "%s"%args.rxnstsv), stdout = subprocess.PIPE)
    # l = sbp.communicate()[0]
    mf = open('/tmp/matched.log','w')
    mmf = open('/tmp/multimatched.log','w')
    nmf = open('/tmp/notmatched.log','w')
    wells = []
    media = collections.defaultdict(list)
    ncpds = {}
    with open(args.mediafile,'r') as infile:
        for line in infile:
            l = line.split('\t')
            n = l[0]
            if args.compoundname not in n:
                continue
            n = n.replace(args.compoundname+'-', '')
            c = l[5:-1]
            sbp = subprocess.Popen(("grep", "-i", "%s"%n, "%s"%args.biolog), stdout = subprocess.PIPE)
            match = (sbp.communicate()[0]).decode('ascii')
            #print(len(match))
            if not match:
                print('?', n, ' was not matched')
                nmf.write(n+'\n')
            else:
                mm = match.split('\n')
                # if len(mm)>2:
                #     mmf.write(n+'\n')
                #     for m in range(len(mm)-1):
                #         x = mm[m].split('\t')
                #         mmf.write('>'+x[0]+'\t'+n+'\t'+x[1]+'\n')
                #     continue
                for m in range(len(mm)-1):
                    x = mm[m].split('\t')
                    if n != x[1]:
                        print('!', n, ' was also matched with ', x[1], x[0], '; this will be skipped')
                        mmf.write('ls '+args.outpath+x[0]+'_'+x[1]+'.tsv\n')
                        continue
                    with open(args.outpath+x[0]+'_'+n+'.tsv','w') as outfile:
                        outfile.write('ID\tminFlux\tmaxFlux\tconcentration\n')
                        for cpd in c:
                            outfile.write(cpd.replace(';','\t')+'\n')
                            media[x[0]].append(cpd.split(';')[0])
                    mf.write(x[0]+'\t'+n+'\n')
                    wells.append(x[0])
    #wells.sort()
    mf.close()
    mmf.close()
    nmf.close()
    #print(media[wells[0]])

    nfw=[]
    #subs = map(lambda x: x[x.index('cpd'):x.index('cpd')+8], s)
    #prods = map(lambda x: x[x.index('cpd'):x.index('cpd')+8], p)
    print('X The following metabolites were not found in the DB (missing ID?)')
    with open(args.biolog,'r') as infile:
        for line in infile:
            l = line.split('\t')
            if l[0] == 'Well':
                continue
            if l[0] not in wells:
                print('\t',l[0], l[1])
                nfw.append(l[0])
    aw = nfw+wells
    aw.sort()
    nfw.sort()
    #print('\n'.join(aw))
    #with open(args.outpath+'biomass.m','w') as outfile:
    #    outfile.write('function [s, p] = biomass()\n')
    #    outfile.write('s = {\''+'\', \''.join(subs)+'\'};\n')
    #    outfile.write('p = {\''+'\', \''.join(prods)+'\'};\n')
    #    outfile.write('end')

    baselinemedia = set(media[wells[0]]).intersection(set(media[wells[1]]))
    print('> Baseline media composition:')
    for m in baselinemedia:
        sbp = subprocess.Popen(("grep", "-i", "%s"%m, "%s"%args.db), stdout = subprocess.PIPE)
        match = (sbp.communicate()[0]).decode('ascii').split('\t')[2]
        print('\t', m, match)
    outfile = open(args.outpath+'ncompounds.tsv','w')
    outfile.write('Well\tSEEDID\n')
    for w in wells:
        ncpds[w] = list(set(media[w]) - baselinemedia)[0]
        outfile.write(w+'\t'+ncpds[w]+'\n')
    outfile.close()
    blm = list(baselinemedia)
    blm.sort()
    outfile = open(args.outpath+'baselinemedia.m','w')
    outfile.write('function [m] = baselinemedia()\n')
    outfile.write('m = {\''+'\', \''.join(blm)+'\'};\n')
    outfile.write('end')
    outfile.close()


def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-m', '--mediafile', help='tsv file with media definition', default='../data/ModelSEEDdata/media_list_with_meta.tsv')
    parser.add_argument('-c', '--compoundname', help='name of source compound', default='Nitrogen')
    parser.add_argument('-b', '--biolog', help='name of biolog tsv file', default='../rhizobiumRoot491/data/biolog_summary.tsv')
    parser.add_argument('-o', '--outpath', help='path for output file', default='../rhizobiumRoot491/media/singleNMedia/')
    parser.add_argument('-d', '--db', help='database file', default="../data/ModelSEEDdata/compounds.tsv")
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
