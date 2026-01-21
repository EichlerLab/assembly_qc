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

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.collections import PatchCollection
import pandas as pd


# args = pd.Series(['chm13', 'h1_merge.bed', 'h1.test.pdf', None, 'h1_test', 'h1_merge.bed'], index=['ref_name', 'asm1', 'outfile', 'chrom', 'sample', 'asm2'])


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for each assembly).  The rest of this script will be prepping data
# for input to this function
#

# --- plotting params ---
# 0: upper / 1: lower
asm_color_dict = {'asm1' : {0 : 'tab:orange', 1 : 'tab:blue'},
                  'asm2' : {0 : 'tab:purple', 1 : 'tab:green'}}

base_chrom_height = 1.25

wide = snakemake.params.wide
sample = snakemake.wildcards.sample
ref = snakemake.params.ref

if wide:
    figsize = (12, 8) # expand width when split == True
    chrom_height = 2
    gene_height  = 3
    gene_padding = 0.5

else:
    figsize = (6,8)
    chrom_height = 1.25  
    gene_height  = 2     
    gene_padding = 0.5

# Keep spacing based on fixed baseline sizes so visible thickness changes
# are not canceled out by proportionally larger spacing.
BASE_CHROM_HEIGHT = 1.25   # baseline for spacing (do not change with args)
BASE_GENE_HEIGHT  = 2.0    # baseline for spacing (do not change with args)

# spacing should NOT scale with gene_height/chrom_height, otherwise thickness changes will look the same proportionally.
chrom_spacing = (BASE_CHROM_HEIGHT + 4*BASE_GENE_HEIGHT + gene_padding) * 1.25

# --- chromosomes to draw ---
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
ideo['width'] = ideo['end'] - ideo['start']


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
ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])
XMAX = float(ideo.groupby('#chrom')['end'].max().max())

# --- layout bookkeeping ---
ybase = 0
chrom_ybase, gene_ybase, chrom_centers = {}, {}, {}
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom]   = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom]    = ybase + chrom_height + gene_padding
    ybase += chrom_height + chrom_spacing

# --- load assemblies, assign vert (0/1 alternating per contig) ---
def load_asm(path, start_vert=0):
    df = pd.read_csv(
        path,
        sep='\t',
        header=None,
        engine='python',
        on_bad_lines='skip',
        usecols=[0,1,2,3,4],
        names=['#chrom','start','end','name','mapq'],
    )
    df = df[df['#chrom'].apply(lambda x: x in chromosome_list)].copy()
    df['width'] = df['end'] - df['start']
    df = df.sort_values(['#chrom','start','end']).copy()

    acc = []
    for chrom in df['#chrom'].unique():
        chrom_df = df.loc[df['#chrom']==chrom].copy()
        contig_dict = {}
        vert = None
        for contig in chrom_df['name']:
            if vert is None:
                vert = start_vert # first contig lane
            else:
                vert = 1 - vert # alternate 0/1
            contig_dict[contig] = vert
        chrom_df['vert'] = chrom_df['name'].map(contig_dict)
        acc.append(chrom_df)
    return pd.concat(acc, ignore_index=True)

asm_dict = {}
# match original staggering: asm1 starts at vert=1, asm2 at vert=0
asm_dict['asm1'] = load_asm(asm1, start_vert=0)
if asm2:
    if wide:
        asm_dict['asm2'] = load_asm(asm2, start_vert=1)
    else:
        asm_color_dict['asm2'] = {1 : 'tab:purple', 0 : 'tab:green'}
        asm_dict['asm2'] = load_asm(asm2, start_vert=0)


# --- drawing helper ---
def chromosome_collections(df, y_positions, height, ptype, pname, *, split=False, keep_vert=None, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes object.
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
    split: boolean
        Only draw rows with df.vert == keep_vert, and place them at a single band (no vertical staggering). Otherwise: stagger by vert (0/1) like original.
            
    Additional kwargs are passed to BrokenBarHCollection    
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']

    for chrom, group in df.groupby('#chrom'):
        xranges = group[['start','width']].values
        pivots  = group[['vert']].values if 'vert' in group else None
        ends    = group[['end']].values
        starts  = group[['start']].values
        names   = group[['name']].values

        if ptype == 'ideo':
            yield BrokenBarHCollection(
                xranges, (y_positions[chrom], height),
                facecolors=group['colors'], **kwargs
            )
            continue

        # aln/fill for contigs
        current_sequence = None
        start_hatch = None
        for i, xval in enumerate(xranges):
            v = int(pivots[i][0]) if pivots is not None else 0
            if split and keep_vert is not None and v != keep_vert:
                # skip the other hap/vert
                continue
            # y-offset control: always stagger by vert (0/1) within each panel
            local_vert_for_offset = v

            if ptype == 'aln':
                if split:
                    # In split mode, draw both asm1 and asm2 above the ideogram (same vertical lane)
                    y0 = y_positions[chrom] + (height * (0)**(local_vert_for_offset)+0.2)
                else:
                    # In non-split mode, keep original layout: asm1 above, asm2 below
                    if pname == 'asm1':
                        y0 = y_positions[chrom] + (height * (0)**(local_vert_for_offset)) + 0.2
                    else:
                        y0 = y_positions[chrom] + (height * (0)**(local_vert_for_offset)) - (gene_height*2 + chrom_height + gene_padding*2) - 0.2

                yield BrokenBarHCollection([xval], (y0, height),
                                           facecolors=asm_color_dict[pname][v], **kwargs)
            else:
                # fill hatchers between adjacent segments of same contig "name"
                sequence = names[i][0]
                if sequence != current_sequence:
                    current_sequence = sequence
                    start_hatch = ends[i][0]
                else:
                    hatch_coords = [start_hatch, starts[i][0] - start_hatch]

                    if pname == 'asm1': # always upper
                        y0 = y_positions[chrom] + (height * (0)**(local_vert_for_offset)) + 0.2
                    else: # asm2
                        if split: # upper
                            y0 = y_positions[chrom] + (height * (0)**(local_vert_for_offset)) + 0.2
                        else: # lower
                            y0 = y_positions[chrom] + (height * (0)**(local_vert_for_offset)) - (gene_height*2 + chrom_height + gene_padding*2) - 0.2
                        
                    yield BrokenBarHCollection([hatch_coords], (y0, height),
                                               facecolors=asm_color_dict[pname][v],
                                               hatch='/', **kwargs)
                current_sequence = sequence
                start_hatch = ends[i][0]

    if del_width:
        del df['width']

# --- figure/axes ---
if wide:
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=figsize, sharey=True, constrained_layout=True)
    axes = [("hap1 (vert=1)", ax_left, 0), ("hap2 (vert=0)", ax_right, 1)]
else:
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    ax = fig.add_subplot(111)

# --- draw frames and ideograms ---
def draw_chrom_frames(ax_):
    for chrom in chromosome_list:
        xmax = float(ideo.loc[ideo['#chrom']==chrom, 'end'].max())
        ax_.barh(y=chrom_ybase[chrom] + (chrom_height/2),
                 width=xmax, height=chrom_height,
                 edgecolor='black', color="None")

def add_ideograms(ax_):
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, 'ideo', 'ideogram', edgecolors="black",linewidths=0.3):
        ax_.add_collection(collection)

if wide:
    for _, ax_, _keep in axes:
        draw_chrom_frames(ax_)
        add_ideograms(ax_)
else:
    draw_chrom_frames(ax)
    add_ideograms(ax)

# --- draw assemblies (contigs) ---
def draw_assemblies(ax_, keep_vert):
    for asm in asm_dict:
        for modifier in ['aln', 'fill']:

            if wide and 'asm2' in asm_dict:
                if asm == 'asm1':
                    for col in chromosome_collections(asm_dict['asm1'], gene_ybase, gene_height, modifier, 'asm1',
                                                    alpha=1.0, linewidths=0.5, split=True):
                        ax_left.add_collection(col)
                elif asm == 'asm2':
                    for col in chromosome_collections(asm_dict['asm2'], gene_ybase, gene_height, modifier, 'asm2',
                                                    alpha=1.0, linewidths=0.5, split=True):
                        ax_right.add_collection(col)
            elif wide and (not 'asm2' in asm_dict):
                for col in chromosome_collections(asm_dict[asm], gene_ybase, gene_height, modifier, asm,
                                                alpha=1.0, linewidths=0.5, split=True, keep_vert=0):
                    ax_left.add_collection(col)
                for col in chromosome_collections(asm_dict[asm], gene_ybase, gene_height, modifier, asm,
                                                alpha=1.0, linewidths=0.5, split=True, keep_vert=1):
                    ax_right.add_collection(col)
            else:
                for col in chromosome_collections(asm_dict[asm], gene_ybase, gene_height, modifier, asm,
                                                alpha=1.0, linewidths=0.5):
                    ax.add_collection(col)


if wide:
    # left: vert==0, right: vert==1
    draw_assemblies(ax_left, 0)
    draw_assemblies(ax_right, 1)
else:
    # as orignal way (staggered by vert)
    draw_assemblies(ax, None)

# --- legend (single) ---
class MulticolorPatch(object):
    def __init__(self, colors):
        self.colors = colors

class MulticolorPatchHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        width, height = handlebox.width, handlebox.height
        patches = []
        for i, c in enumerate(orig_handle.colors):
            patches.append(plt.Rectangle([width/len(orig_handle.colors)*i - handlebox.xdescent,
                                          -handlebox.ydescent],
                                         width/len(orig_handle.colors),
                                         height,
                                         facecolor=c,
                                         edgecolor='none'))
        patch = PatchCollection(patches, match_original=True)
        handlebox.add_artist(patch)
        return patch

if wide:
    target_ax = ax_right 
else:
    target_ax = ax

h, l = target_ax.get_legend_handles_labels()
h.append(MulticolorPatch(['tab:orange', 'tab:blue']));  l.append("asm1")
if 'asm2' in asm_dict:
    h.append(MulticolorPatch(['tab:green', 'tab:purple'])); l.append("asm2")

leg = target_ax.legend(h, l, loc='lower right', fontsize=10,
                       labelspacing=0.6, borderpad=0.5, handlelength=2.0,
                       handler_map={MulticolorPatch: MulticolorPatchHandler()}
                    )

# --- ticks / title / save ---
def finalize_axes(ax_):
    ax_.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax_.set_yticklabels(chromosome_list)

if wide:
    finalize_axes(ax_left)
    finalize_axes(ax_right)
    ax_left.margins(x=0.0)
    ax_right.margins(x=0.0)
    ax_right.tick_params(axis='y', which='both', labelleft=True, left=True)
    ax_right.spines['left'].set_visible(True)
    ax_left.set_xlim(0, XMAX * 1.05)   # small pad avoids right-edge clipping
    ax_right.set_xlim(0, XMAX * 1.05)    
    ax_left.set_title("Haplotype 1", fontsize=10)
    ax_right.set_title("Haplotype 2", fontsize=10)
    fig.suptitle(f"{sample} → {ref}", fontsize=15)
    fig.set_constrained_layout_pads(w_pad=0.0, h_pad=0.0, wspace=0.02, hspace=0.0)
else:
    finalize_axes(ax)
    ax.set_title(f"{sample} → {ref}")

plt.savefig(snakemake.output.pdf, bbox_inches='tight', dpi=600)
plt.close('all')
