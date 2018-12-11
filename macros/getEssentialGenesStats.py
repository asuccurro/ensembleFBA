#!/usr/bin/python3
#************************************
#**  author: Antonella Succurro    **
#**  email:a.succurro[AT]gmail.com **
#**                                **
#**  created:       2018/05/24     **
#**  last modified: 2018/05/24     **
#************************************

import argparse
import pandas
import numpy as np
import string
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from makeBiologFigure import getNetworkStats, addMajorityCol
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import venn
import itertools

mylabels = {'ammo': 'Ammonia', 'lysi': 'L-Lysine', 'seri': 'L-Serine', 'urea': 'Urea', 'glut': 'L-Glutamic-Acid',
            'Root9': 'Root9', 'Root491':'Root491', 'Root66D1': 'Root66D1'}
mycolors={'lysi': '#8DA0CB',
          'glut': '#FC8D62',
          'seri': '#E78AC3',
          'urea': '#a2cb8d',
          'ammo': '#66C2A5',
          'Root9': '#e14821',
          'Root491':'#9c21e1',
          'Root66D1': '#dfe121'}

patricids = {'Root9': 'fig|1736604.3', 'Root491':'fig|1736548.3', 'Root66D1': 'fig|1736582.3'}
genbankids= {'Root9': 'ASE33_', 'Root491':'ASD46_', 'Root66D1': 'ASE09_'}


def main():
    args = options()
    verbose = args.verbose

    Nthr = -1
    Nstr = ''
    
    if args.any:
        Nthr = 0
        Nstr = '_any'
    elif args.una:
        Nthr = 50
        Nstr = '_all'
    elif args.maj:
        Nthr = 25
        Nstr = '_maj'
    else:
        print('No threshold chosen, exiting. Please select "any" -A; "unanimity" -U; "majority" -M.')
        return

    sets_ge = {}
    dfs_ge = {}
    
    orgs = 'Root9 Root491 Root66D1'.split(' ')
    cond = 'ammo glut lysi seri urea'.split(' ')

    cmap_c = LinearSegmentedColormap.from_list('cmap_c', [mycolors[x] for x in cond], N=len(cond))
    cmap_o = LinearSegmentedColormap.from_list('cmap_o', [mycolors[x] for x in orgs], N=len(cond))
        
    if args.cond:
        for o in orgs:
            sets_ge[o] = {}
            dfs_ge[o] = {}
            for c in cond:
                #print('%s%s_%s.csv' % (args.iopath, o, c))
                dfs_ge[o][c] = pandas.read_csv('%s%s_%s.csv' % (args.iopath, o, c)).rename(columns={"Row": "genes", "N": c}).set_index('genes')
                dfs_ge[o][c] = dfs_ge[o][c][dfs_ge[o][c][c] > Nthr]
                sets_ge[o][c] = set((dfs_ge[o][c].index))
            if args.distr:
                df = pandas.concat([dfs_ge[o][cond[0]], dfs_ge[o][cond[1]], dfs_ge[o][cond[2]], dfs_ge[o][cond[3]], dfs_ge[o][cond[4]] ], axis=1, sort=True)
                #print(df.head())
                plotBarh(df, cmap_c, '%sessential_genes_distr_%s%s.png' % (args.iopath, o, Nstr), o)
                #plotBar(df, cond, '/tmp/plot_%s.png' % o)
                #stacked_bar_chart(df, cond, 'genes', 'xyz', 'x', 'y', '/tmp/plot_%s.png' % o, Color('#dfe121'), Color('#e14821'))
            getUnique(dfs_ge[o], cond, '%sunique_eg_%s%s.csv' % (args.iopath, o, Nstr), genbankids[o])
            #checkSets(dfs_ge[o], cond, '%sshared_eg_%s%s_' % (args.iopath, o, Nstr))
        plotVenns(sets_ge, orgs, cmap_c, '%svenn_strains%s.png' % (args.iopath, Nstr))
    
    elif args.orgs:
        for c in cond:
            sets_ge[c] = {}
            dfs_ge[c] = {}
            for o in orgs:
                #print('%s%s_%s.csv' % (args.iopath, o, c))
                dfs_ge[c][o] = pandas.read_csv('%s%s_%s.csv' % (args.iopath, o, c)).rename(columns={"Row": "genes", "N": o}).set_index('genes')
                dfs_ge[c][o] = dfs_ge[c][o][dfs_ge[c][o][o] > Nthr]
                sets_ge[c][o] = set((dfs_ge[c][o].index))
            if args.distr:
                df = pandas.concat([dfs_ge[c][orgs[0]], dfs_ge[c][orgs[1]], dfs_ge[c][orgs[2]] ], axis=1, sort=True)
                #print(df.head())
                plotBarh(df, cmap_o, '%splot_%s%s.png' % (args.iopath, c, Nstr))
            checkSets(dfs_ge[c], orgs, '%sshared_eg_%s%s_' % (args.iopath, c, Nstr))
        plotVenns(sets_ge, cond, cmap_o, '%svenn_conditions%s.png' % (args.iopath, Nstr))


    return dfs_ge

def checkSets(df, ls, oname):
    allpairs = list(itertools.combinations(ls, 2))
    for p in allpairs:
        idx = set(df[p[0]].index) & set(df[p[1]].index)
        tmp = pandas.concat([df[p[0]].loc[list(idx)], df[p[1]].loc[list(idx)]], axis=1, sort=True)
        tmp.to_csv('%s%s.csv' % (oname, '_'.join(p)))
    return


def getUnique(df, ls, oname, strainpre):
    for i in range(len(ls)):
        allother = ls[:]
        cond = allother.pop(i)
        sao = set()
        print(cond)
        for j in allother:
            sao = sao | set(df[j].index)
        uniques = set(df[cond].index) - sao
        with open(oname, 'w') as ofile:
            for u in uniques:
                ofile.write('%s%s\n' % (strainpre, u))
        
    return

def plotVenns(df, subs, mycmap, oname):

    fig, axes = plt.subplots(len(subs), 1, figsize=(5, 4*len(subs)))
    for i, c in enumerate(subs):
        axes[i].set_title(mylabels[c])
        venn.venn(df[c], cmap=mycmap, ax=axes[i], fontsize=6, legend_loc="best")
    fig.tight_layout(pad=0.1)
    plt.savefig(oname)
    return

def plotVenns4(df, subs, mycmap, oname):

    rw = [0, 0, 1, 1]
    cl = [0, 1, 0, 1]
    fig, axes = plt.subplots(2, 2)
    for i, c in enumerate(subs):
        axes[rw[i],cl[i]].set_title(mylabels[c])
        venn.venn(df[c], cmap=mycmap, ax=axes[rw[i],cl[i]], fontsize=6, legend_loc="best")
    fig.tight_layout(pad=0.1)
    plt.savefig(oname)
    return

def plotBarh(df, mycmap, oname, strain):
    f, ax = plt.subplots(figsize=(6, 15))
    df['total'] = df.sum(axis=1)
    df=df.sort_values("total", ascending=False)
    df = df.drop('total', axis=1)
    df.plot.barh(stacked=True, colormap=mycmap, ax=ax)
    ax.set_title('Predicted essential genes for '+strain)
    ax.set_xlabel('N of networks')
    plt.savefig(oname)

def plotBar(df, subbars, oname):

    sns.set(style="whitegrid")

    f, ax = plt.subplots(figsize=(6, 15))

    df['total'] = df.sum(axis=1)
    df=df.sort_values("total", ascending=False)
    df['genes'] = df.index
    
    sns.set_color_codes("pastel")
    #sns.barplot(x="total", y="genes", data=df, label="Total", color="b")

    for sb in subbars:
        print(sb, mylabels[sb], mycolors[sb])
        sns.barplot(x=sb, y="genes", data=df, label=mylabels[sb], color=mycolors[sb])

    ax.legend(ncol=len(subbars), loc="lower right", frameon=True)
    ax.set(xlim=(0, 24), ylabel="yyy", xlabel="bbb")
    sns.despine(left=True, bottom=True)
    plt.savefig(oname)

def stacked_bar_chart(pivoted_df, stack_vals, level_values_field, chart_title, x_label, y_label, filename, color1, color2):
    # https://gist.github.com/extrospective/0f4fe69304184d813f982035d9684452
    # stacked_bar_chart: draws and saves a barchart figure to filename
    #
    # pivoted_df: dataframe which has been pivoted so columns correspond to the values to be plotted
    # stack_vals: the column names in pivoted_df to plot
    # level_values_field: column in the dataframe which has the values to be plotted along the x axis (typically time dimension)
    # chart_title: how to title chart
    # x_label: label for x axis
    # y_label: label for y axis
    # filename: full path filename to save file
    # color1: first color in spectrum for stacked bars
    # color2: last color in spectrum for stacked bars; routine will select colors from color1 to color2 evenly spaced
    #
    # Implementation: based on (http://randyzwitch.com/creating-stacked-bar-chart-seaborn/; https://gist.github.com/randyzwitch/b71d47e0d380a1a6bef9)
    # this routine draws overlapping rectangles, starting with a full bar reaching the highest point (sum of all values), and then the next shorter bar
    # and so on until the last bar is drawn.  These are drawn largest to smallest with overlap so the visual effect is that the last drawn bar is the
    # bottom of the stack and in effect the smallest rectangle drawn.
    #
    # Here "largest" and "smallest" refer to relationship to foreground, with largest in the back (and tallest) and smallest in front (and shortest).
    # This says nothing about which part of the bar appear large or small after overlap.
    #
    color_spectrum = list(color1.range_to(color2, len(stack_vals)))
    plt.clf()
    #
    stack_total_column = 'Stack_subtotal_xyz'  # placeholder name which should not exist in pivoted_df
    bar_num = 0
    legend_rectangles = []
    legend_names = []
    for bar_part in stack_vals:    # for every item in the stack we need to compute a rectangle
        #stack_color = color_spectrum[bar_num].get_hex_l()  # get_hex_l ensures full hex code of color
        stack_color = mycolors[bar_part]
        sub_count = 0
        pivoted_df[stack_total_column] = 0
        stack_value = ""
        for stack_value in stack_vals:  # for every item in the stack we create a new subset [stack_total_column] of 1 to N of the sub values
            pivoted_df[stack_total_column] += pivoted_df[stack_value]  # sum up total
            sub_count += 1
            if sub_count >= len(stack_vals) - bar_num:  # we skip out after a certain number of stack values
                break
        # now we have set the subtotal and can plot the bar.  reminder: each bar is overalpped by smaller subsequent bars starting from y=0 axis
        bar_plot = sns.barplot(data=pivoted_df, x=pivoted_df.index.get_level_values(level_values_field),
                               y=stack_total_column, color=stack_color)
        legend_rectangles.append(plt.Rectangle((0,0),1,1,fc=stack_color, edgecolor = 'none'))
        legend_names.append(stack_value)   # the "last" stack_value is the name of that part of the stack
        bar_num += 1
    l = plt.legend(legend_rectangles, legend_names, loc=2, ncol = 1, prop={'size':12})
    l.draw_frame(False)
    bar_plot.set(xlabel=x_label, ylabel=y_label)
    plt.tight_layout()
    plt.title(chart_title)
    sns.despine(left=True)
    plt.savefig(filename)
    
    
def options():
    '''define here in-line arguments'''
    parser = argparse.ArgumentParser(description='Parsing options')
    parser.add_argument('-V', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-D', '--distr', help='plot distribution', action='store_true')
    parser.add_argument('-A', '--any', help='use any threshold', action='store_true')
    parser.add_argument('-U', '--una', help='use unanimity threshold', action='store_true')
    parser.add_argument('-M', '--maj', help='use majority threshold', action='store_true')
    parser.add_argument('-O', '--orgs', help='compare organisms across conditions', action='store_true')
    parser.add_argument('-C', '--cond', help='compare conditions across organisms', action='store_true')
    parser.add_argument('-p', '--iopath', help='path for input and output file', default='../outputs/geneEssentiality/')
    args = parser.parse_args()
    if args.verbose:
        print("verbosity turned on")
        print(args)
    return args

if __name__=="__main__":
    main()

