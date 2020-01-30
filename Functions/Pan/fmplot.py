#!/usr/bin/env python

__version__ = '0.1.0'

def get_options():
    import argparse

    # create the top-level parser
    description = "Create plots from roary outputs"
    parser = argparse.ArgumentParser(description = description,
                                     prog = 'fmplot.py')

    parser.add_argument('tree', action='store',
                        help='Newick Tree file', default='accessory_binary_genes.fa.newick')
    parser.add_argument('spreadsheet', action='store',
                        help='Roary gene presence/absence spreadsheet', default='gene_presence_absence.csv')

    parser.add_argument('--labels', action='store_true',
                        default=False,
                        help='Add node labels to the tree (up to 18 chars)')
    parser.add_argument('--format',
                        choices=('png',
                                 'tiff',
                                 'pdf',
                                 'svg'),
                        default='pdf',
                        help='Output format [Default: pdf]')
    parser.add_argument('-N', '--skipped-columns', action='store',
                        type=int,
                        default=14,
                        help='First N columns of Roary\'s output to exclude [Default: 14]')
    
    parser.add_argument('--version', action='version',
                         version='%(prog)s '+__version__)

    return parser.parse_args()

if __name__ == "__main__":
    options = get_options()

    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style('white')

    import os
    import pandas as pd
    import numpy as np
    from Bio import Phylo

    t = Phylo.read(options.tree, 'newick')

    # Max distance to create better plots
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    # Load roary
    roary = pd.read_csv(options.spreadsheet,
                         sep=',',
                         low_memory=False)
    # Set index (group name)
    roary.set_index('Gene', inplace=True)
    # Drop the other info columns
    roary.drop(list(roary.columns[:options.skipped_columns-1]), axis=1, inplace=True)

    # Transform it in a presence/absence matrix (1/0)
    roary.replace('.{2,100}', 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)

    # Sort the matrix by the sum of strains presence
    idx = roary.sum(axis=1).sort_values(ascending=False).index
    roary_sorted = roary.loc[idx]

    # Pangenome frequency plot
    plt.figure(figsize=(7, 5))
    #histtype: bar, barstacked, step, stepfilled
    plt.hist(roary.sum(axis=1), roary.shape[1],
             histtype="stepfilled", alpha=.7, facecolor="Teal", edgecolor="black")

    plt.xlabel('No. of genomes')
    plt.ylabel('No. of genes')

    sns.despine(left=True,
                bottom=True)
    plt.savefig('pangenome_frequency.%s'%options.format, dpi=300)
    plt.clf()

    # Sort the matrix according to tip labels in the tree
    roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]

    # Plot presence/absence matrix against the tree
    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))

        ax1=plt.subplot2grid((1,40), (0, 15), colspan=30)
        a=ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues,
                   vmin=0, vmax=1,
                   aspect='auto',
                   interpolation='none',
                    )
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.axis('off')

        ax = fig.add_subplot(1,2,1)
        # matplotlib v1/2 workaround
        try:
            ax=plt.subplot2grid((1,40), (0, 0), colspan=10, facecolor='white')
        except AttributeError:
            ax=plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('Gene matrix\n(%d gene clusters)'%roary.shape[0])

        if options.labels:
            fsize = 12 - 0.1*roary.shape[1]
            if fsize < 5:
                fsize = 5
            with plt.rc_context({'font.size': fsize}):
                Phylo.draw(t, axes=ax, 
                           show_confidence=False,
                           label_func=lambda x: str(x)[:18],
                           xticks=([],), yticks=([],),
                           ylabel=('',), xlabel=('',),
                           xlim=(-mdist*0.1,mdist+mdist*0.45-mdist*roary.shape[1]*0.001),
                           axis=('off',),
                           title=('Tree\n(%d strains)'%roary.shape[1],), 
                           do_show=False,
                          )
        else:
            Phylo.draw(t, axes=ax, 
                       show_confidence=False,
                       label_func=lambda x: None,
                       xticks=([],), yticks=([],),
                       ylabel=('',), xlabel=('',),
                       xlim=(-mdist*0.1,mdist+mdist*0.1),
                       axis=('off',),
                       title=('Tree\n(%d strains)'%roary.shape[1],),
                       do_show=False,
                      )
        plt.savefig('pangenome_matrix.%s'%options.format, dpi=300)
        plt.clf()
