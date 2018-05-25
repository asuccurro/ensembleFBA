#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/25     **
#**  last modified: 2018/05/25     **
#************************************

import argparse
import subprocess

def main():
    args = options()
    verbose = args.verbose
    # sbp = subprocess.Popen(("grep", "-w", "%s"%args.bmname, "%s"%args.rxnstsv), stdout = subprocess.PIPE)
    # l = sbp.communicate()[0]
    with open(args.mediafile,'r') as infile:
        for line in infile:
            l = line.split('\t')
            n = l[0]
            if args.compoundname not in n:
                continue
            n = n.replace(args.compoundname+'-', '')
            c = l[5:-1]
            sbp = subprocess.Popen(("grep", "-i", "%s"%n, "../rhizobiumRoot491/data/biolog_summary.csv"), stdout = subprocess.PIPE)
            match = (sbp.communicate()[0]).decode('ascii')
            #print(len(match))
            if not match:
                print('***** \t', n, ' not found')
            else:
                mm = match.split('\n')
                if len(mm)>2:
                    print('***** \t', n, ' has more matches!')
                    for m in range(len(mm)-1):
                        x = mm[m].split('\t')
                        print('***** \t \t', n, x[1], x[0])
                    continue
                for m in range(len(mm)-1):
                    x = mm[m].split('\t')
                    print(n, x[1], x[0])
                    with open(args.outpath+x[0]+'_'+n+'.tsv','w') as outfile:
                        outfile.write('ID\tminFlux\tmaxFlux\tconcentration\n')
                        for cpd in c:
                            outfile.write(cpd.replace(';','\t')+'\n')

                    #print(n, x)


    #subs = map(lambda x: x[x.index('cpd'):x.index('cpd')+8], s)
    #prods = map(lambda x: x[x.index('cpd'):x.index('cpd')+8], p)

    #with open(args.outpath+'biomass.m','w') as outfile:
    #    outfile.write('function [s, p] = biomass()\n')
    #    outfile.write('s = {\''+'\', \''.join(subs)+'\'};\n')
    #    outfile.write('p = {\''+'\', \''.join(prods)+'\'};\n')
    #    outfile.write('end')


def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-m', '--mediafile', help='tsv file with media definition', default='../data/ModelSEEDdata/media_list_with_meta.tsv')
    parser.add_argument('-c', '--compoundname', help='name of source compound', default='Nitrogen')
    parser.add_argument('-o', '--outpath', help='path for output file', default='../rhizobiumRoot491/data/singleNMedia/')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
