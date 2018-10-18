#!../venvpy/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/29     **
#**  last modified: 2018/05/29     **
#************************************

import argparse
import pandas

def main():
    args = options()
    verbose = args.verbose
    fname = '../data/MPI%s/MPI%s-reactions.tsv' % (args.id, args.id)
    df = pandas.read_csv(fname, sep='\t')
    df.to_csv(fname, sep='\t', index=False)

def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-i', '--id', help='id for the organism', default='Root9')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()
