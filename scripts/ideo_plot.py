#!/bin/env python


"""
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.
1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::
    #chrom  chromStart  chromEnd  name    gieStain
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg
2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.
"""
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.collections import PatchCollection
import pandas as pd
import numpy as np

# args = pd.Series(['chm13', 'h1_merge.bed', 'h1.test.pdf', None, 'h1_test', 'h1_merge.bed'], index=['ref_name', 'asm1', 'outfile', 'chrom', 'sample', 'asm2'])


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for each assembly).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height, ptype, pname, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    # plot each df with 
    for chrom, group in df.groupby('#chrom'):
        xranges = group[['start', 'width']].values

        # if xranges.ndim == 1 or xranges.shape[1] != 2:
        #     xranges = np.array([x if isinstance(x, (list, np.ndarray)) and len(x) == 2 else [0, 0] for x in xranges])
        if ptype == 'aln':
            pivot = group[['vert']].values
            for i, xval in enumerate(xranges):

                if pname == 'asm1':
                    yield BrokenBarHCollection(
                        [xval], (y_positions[chrom]+(height*(0)**(pivot[i][0])), height), facecolors=asm_color_dict[pname][pivot[i][0]], **kwargs
                    )
                else:
                    yield BrokenBarHCollection(
                        [xval], (y_positions[chrom]+(height*(0)**(pivot[i][0]))-(gene_height*2+chrom_height+gene_padding*2), height), facecolors=asm_color_dict[pname][pivot[i][0]], **kwargs
                    )
        elif ptype == 'ideo':
            yield BrokenBarHCollection(
                    xranges, (y_positions[chrom], height), facecolors=group['colors'], **kwargs
                )
        else:
            ends = group[['end']].values
            starts = group[['start']].values
            pivot = group[['vert']].values
            contigs = group[['name']].values
            current_sequence = None
            for i, xval in enumerate(xranges):
                sequence = contigs[i]
                if sequence != current_sequence:
                    current_sequence = sequence
                    start_hatch = ends[i]
                elif current_sequence is not None:
                    # ax.barh(current_sequence, start - start_hatch, left=start_hatch, height=0.5,
                            # color=sequence_colors[current_sequence], hatch='/')
                    hatch_coords = [start_hatch, starts[i]-start_hatch]
                    if pname == 'asm1':
                        yield BrokenBarHCollection(
                            [hatch_coords], (y_positions[chrom]+(height*(0)**(pivot[i][0])), height), facecolors=asm_color_dict[pname][pivot[i][0]], hatch='/', **kwargs
                        )
                    else:
                        yield BrokenBarHCollection(
                            [hatch_coords], (y_positions[chrom]+(height*(0)**(pivot[i][0]))-(gene_height*2+chrom_height+gene_padding*2), height), facecolors=asm_color_dict[pname][pivot[i][0]], hatch='/', **kwargs
                        )
                current_sequence = sequence
                start_hatch = ends[i]
    if del_width:
        del df['width']


asm_color_dict = {'asm1' : {0 : 'tab:blue', 1 : 'tab:orange'}, 'asm2' : {0 : 'tab:green', 1 : 'tab:purple'}}

# Height of each ideogram
chrom_height = 1.25

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 2

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.5


# Spacing between consecutive ideograms
chrom_spacing = (chrom_height+gene_height*4+gene_padding)*1.25
# chrom_spacing = 8

# Width, height (in inches)
figsize = (6, 8)



# Decide which chromosomes to use
with open(snakemake.params.chroms) as finp:
    token = finp.read().strip().split("\n")
chromosome_list = [ line.split()[0] for line in token if not line.startswith("#") ]

all_beds = snakemake.input
asm1 = all_beds[0]
if len(all_beds) == 2:
    asm2 = all_beds[1]
else:
    asm2 = False

cytoband = snakemake.params.cyto
if not cytoband == "Null":
    ideo = pd.read_csv(
        cytoband,
        sep='\t'
    )
else:
    fai = snakemake.params.path+".fai"
    ideo = pd.read_csv(fai, sep="\t", header=None, usecols=[0,1], names=["#chrom","end"])
    ideo["start"] = 0
    ideo["name"] = ideo["#chrom"]
    ideo["gieStain"] = "gpos25"


ideo = ideo[ideo['#chrom'].apply(lambda x: x in chromosome_list)]
ideo["start"] = ideo["start"].astype("int")
ideo["end"] = ideo["end"].astype("int")

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase + chrom_height + gene_padding
    ybase += chrom_height + chrom_spacing


# Add a new column for width
ideo['width'] = ideo['end'] - ideo['start']

# Colors for different chromosome stains
color_lookup = {
    'gneg': (1., 1., 1.),
    'gpos25': (.6, .6, .6),
    'gpos50': (.4, .4, .4),
    'gpos75': (.2, .2, .2),
    'gpos100': (0., 0., 0.),
    'acen': (.8, .4, .4),
    'gvar': (.8, .8, .8),
    'stalk': (.9, .9, .9),
}

# Add a new column for colors
ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])


# Same thing for genes
asm_dict = {}

df = pd.read_csv(asm1, names=['#chrom', 'start', 'end', 'name', 'mapq'], sep='\t', header=None)
df = df[df['#chrom'].apply(lambda x: x in chromosome_list)]
df['width'] = df['end'] - df['start']
df = df.sort_values(['#chrom', 'start', 'end']).copy()

chrom_all = pd.DataFrame()
for chrom in df['#chrom'].unique():
    chrom_df = df.loc[df['#chrom'] == chrom].copy()
    contig_dict = {}
    vert = None
    for i, contig in enumerate(chrom_df['name']):
        vert = 1 if vert != 1 else 0
        contig_dict[contig] = vert
    chrom_df['vert'] = chrom_df['name'].apply(lambda val: contig_dict[val])
    chrom_all = pd.concat([chrom_all, chrom_df])
asm_dict['asm1'] = chrom_all.copy()

# Load and transform 
if asm2:
    df = pd.read_csv(asm2, names=['#chrom', 'start', 'end', 'name', 'mapq'], sep='\t', header=None)
    df = df[df['#chrom'].apply(lambda x: x in chromosome_list)]
    df['width'] = df['end'] - df['start']
    df = df.sort_values(['#chrom', 'start', 'end']).copy()
    chrom_all = pd.DataFrame()
    for chrom in df['#chrom'].unique():
        chrom_df = df.loc[df['#chrom'] == chrom].copy()
        contig_dict = {}
        vert = None
        for i, contig in enumerate(chrom_df['name']):
            vert = 0 if vert != 0 else 1
            contig_dict[contig] = vert
        chrom_df['vert'] = chrom_df['name'].apply(lambda val: contig_dict[val])
        chrom_all = pd.concat([chrom_all, chrom_df])
    asm_dict['asm2'] = chrom_all.copy()


fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)


for i, chrom in enumerate(chromosome_list):
    ax.barh(y=chrom_ybase[chrom]+(chrom_height/2), width=max(ideo.loc[ideo['#chrom'] == chrom]['end']), height=chrom_height, edgecolor='black', color="None")

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, 'ideo', 'ideogram'):
    ax.add_collection(collection)

# ...and the gene data
for asm in asm_dict:
    print(f"plotting {asm}")
    for modifier in ['aln', 'fill']:

        for collection in chromosome_collections(
            asm_dict[asm], gene_ybase, gene_height, modifier, asm, alpha=1.0, linewidths=0.5,
        ):
            ax.add_collection(collection)

# define an object that will be used by the legend
class MulticolorPatch(object):
    def __init__(self, colors):
        self.colors = colors
        
# define a handler for the MulticolorPatch object
class MulticolorPatchHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        width, height = handlebox.width, handlebox.height
        patches = []
        for i, c in enumerate(orig_handle.colors):
            patches.append(plt.Rectangle([width/len(orig_handle.colors) * i - handlebox.xdescent, 
                                          -handlebox.ydescent],
                           width / len(orig_handle.colors),
                           height, 
                           facecolor=c, 
                           edgecolor='none'))

        patch = PatchCollection(patches,match_original=True)

        handlebox.add_artist(patch)
        return patch


h, l = ax.get_legend_handles_labels()

h.append(MulticolorPatch(['tab:orange', 'tab:blue']))
l.append("asm1")

h.append(MulticolorPatch(['tab:green', 'tab:purple']))
l.append("asm2")

plt.legend(h, l, loc='lower right', fontsize="small",
         handler_map={MulticolorPatch: MulticolorPatchHandler()}) 
         # bbox_to_anchor=(.125,.875))

# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.set_title(f"{snakemake.wildcards.asm} $\\rightarrow$ {snakemake.params.ref}")
#ax.axis('off')
plt.savefig(snakemake.output.pdf, bbox_inches='tight')
plt.close('all')
