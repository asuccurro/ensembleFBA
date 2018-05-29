#!/usr/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/24     **
#**  last modified: 2018/05/24     **
#************************************

import argparse
import subprocess

def main():
    args = options()
    verbose = args.verbose
    sbp = subprocess.Popen(("grep", "-w", "%s"%args.bmname, "%s"%args.rxnstsv), stdout = subprocess.PIPE)
    l = sbp.communicate()[0].decode('ascii')
    #with open('bio1.txt','r') as infile:
    #    for line in infile:
    #        l = line

    e = l.split('\t')[8].split('=>')
    s=e[0].split('+')
    p=e[1].split('+')

    subs = map(lambda x: x[x.index('cpd'):x.index('cpd')+8], s)
    prods = map(lambda x: x[x.index('cpd'):x.index('cpd')+8], p)

    with open(args.outpath+'biomass.m','w') as outfile:
        outfile.write('function [s, p] = biomass()\n')
        outfile.write('s = {\''+'\', \''.join(subs)+'\'};\n')
        outfile.write('p = {\''+'\', \''.join(prods)+'\'};\n')
        outfile.write('end')


def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-r', '--rxnstsv', help='path to the GSMNM reactions tsv file', default='../rhizobiumRoot491/draftgsmn/tsv/rhizobiumRoot491-reactions.tsv')
    parser.add_argument('-b', '--bmname', help='biomass reaction name', default='bio1')
    parser.add_argument('-o', '--outpath', help='path for output file', default='../rhizobiumRoot491/draftgsmn/')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
