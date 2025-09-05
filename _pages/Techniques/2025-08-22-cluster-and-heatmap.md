---
title: "分子动力学聚类分析与热图可视化技术"
date: "2025-08-22"
tags: [clustering, heatmap, molecular-dynamics, data-analysis, visualization, conformational-analysis, trajectory-analysis, vmd, python, data-visualization]
---

# 分子动力学聚类分析与热图可视化技术

## 聚类分析

首先需要准备一个包含每个聚类中帧ID的 `clus_result.dat` 文件，格式如下（每个聚类的第一个数字是中心构象）：

```
cluster 1: 3722
3946 0 1 4 10 23 33 36 41 45 46 47 51 54 59 61 62 63 66 67 69 76 80 84 85 ......

cluster 2: 489
1886 2 3 5 8 9 11 13 14 16 17 18 19 20 21 22 24 25 27 30 31 32 34 35 37 38 39 40 42 43 44 48 49....

....
```

在VMD中通过以下TCL脚本生成：

```tcl
# http://github.com/anjibabuIITK/CLUSTER-ANALYSIS-USING-VMD-TCL
set number 9	;# number of clusters, others are tagged 'other'
set rcutoff 1.5  ;# RMSD cutoff. unit: angstrom
set step_size 1
set nframes [molinfo top get numframes]
set inf 0
set nf $nframes
set totframes [expr $nf - 1 ]
set selA [atomselect top "fragment 1 and resid 149 to 156 and backbone"]    ;# select the ligand
set lists [measure cluster $selA num $number cutoff $rcutoff first $inf last $totframes step $step_size distfunc rmsd weight mass]

set file [open "clus_result.dat" w]
for {set i 1} {$i <= [llength $lists]} {incr i} {
    set lst [lindex $lists [expr $i-1]]
    puts $file [format "cluster %d: %d" $i [llength $lst]]
    puts $file $lst
    puts $file ""
}
close $file

# save the coordinates of centroid structures
set c01 [lindex [lindex $lists 0] 0]
set sel [atomselect top all frame $c01]
set real_frame [expr $c01+1]
$sel writegro c01_${real_frame}.gro
puts [format "write the centroid of 1st cluster: frame %d" $real_frame]

set c02 [lindex [lindex $lists 1] 0]
set sel [atomselect top all frame $c02]
set real_frame [expr $c02+1]
$sel writegro c02_${real_frame}.gro
puts [format "write the centroid of 2nd cluster: frame %d" $real_frame]
```

然后使用Python进行可视化：

```python
import matplotlib.pyplot as plt
import numpy as np
import os


def read_vmd_clus_result(file):
    data = []
    with open(file, 'r') as f:
        while f.readline().strip().startswith('cluster'):
            line = f.readline().strip()
            data.append([int(fr) for fr in line.split()])
            _ = f.readline()  # empty
    return data

def get_id_with_time(data):
    # data: output from read_vmd_clus_result()
    # return: a list of tuples, (frame_id, cluster_id)
    # cluster_id starts from 1
    id_with_time = []
    for i in range(len(data)):
        cl = data[i]
        id_with_time += [(fr, i + 1) for fr in cl]
    id_with_time.sort(key=lambda x: x[0])
    return id_with_time

font_le = {'family': 'Times New Roman', 'weight': 'demibold', 'size': 16}
font_la = {'family': 'Times New Roman', 'fontname': 'Times New Roman', 'weight': 'demibold', 'size': 24}
font_tc = {'family': 'Times New Roman', 'fontname': 'Times New Roman', 'weight': 'demibold', 'size': 20}
font_ti = {'family': 'Times New Roman', 'fontname': 'Times New Roman', 'weight': 'demibold', 'size': 28}
font_hu = {'family': 'Times New Roman', 'fontname': 'Times New Roman', 'weight': 'demibold', 'size': 36}

# a framework of the plot
def plot_common(xlabel, ylabel, thickness=2, title=None, size=(8,6), xpad=6, ypad=0, title_pad=0, ticks_size=16, tight=False, ax_color='black'):
    fig, ax = plt.subplots(figsize=size)  # fix the xlabel overflow problem
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(thickness)
    ax.tick_params(width=thickness)
    ax.tick_params(axis='y', colors=ax_color, labelcolor=ax_color)
    # plt.xticks(font='Arial', size=16, weight='bold')
    # plt.yticks(font='Arial', size=16, weight='bold')
    plt.xticks(font='Times New Roman', size=ticks_size, weight='demibold')
    plt.yticks(font='Times New Roman', size=ticks_size, weight='demibold')
    plt.xlabel(xlabel, fontdict=font_la, labelpad=xpad)
    plt.ylabel(ylabel, fontdict=font_la, labelpad=ypad, color=ax_color)
    if title is not None:
        plt.title(title, fontdict=font_ti, pad=title_pad)
    if tight:
        plt.tight_layout()
    return fig, ax

def plot_clustering_id_with_time(idxs, nsperframe, biggest=10, path=None, point=False, size=(8,6), ssize=1):
    # plot the frame_id with cluster_id. Marking the selected centroid frame (point) with a star.
    # nsperframe: convert frame_id to nanosecond
    # biggest: biggest cluster_id shown. Other frames are tagged 'other'.
    plot_common(xlabel='Time (ns)', ylabel='Cluster ID', size=size)
    biggest = min(biggest, int(max(idxs)))
    plt.yticks(np.arange(biggest+1), labels=np.arange(biggest).tolist()+['Other'])

    x = np.arange(len(idxs))*nsperframe
    y = [min(i, biggest) for i in idxs]
    plt.scatter(x, y, s=ssize)
    if point:
        plt.scatter(point*nsperframe, idxs[point], marker='*', s=50, color='r')
    print("The number of clusters: {0:d}".format(len(set(idxs))))
    print("The biggest cluster lasted for {0:.1f} ns ({1:.1%})".format(np.sum(idxs==1)*nsperframe, np.sum(idxs==1)/len(idxs)))
    print("The second  cluster lasted for {0:.1f} ns ({1:.1%})".format(np.sum(idxs==2)*nsperframe, np.sum(idxs==2)/len(idxs)))
    print("The unclustered frames (>=no. {2:d}) occupies {0:.1f} ns or {1:.1%}".format(np.sum(idxs>=biggest)*nsperframe, np.sum(idxs>=biggest)/len(idxs), biggest))
    print("The centroid of the biggest cluster is at {0:.1f} ns.".format(point*nsperframe))
    if path is not None:
        plt.savefig(os.path.join(path,'cluster.png'))
    plt.show()


path = 'xxxxxxxxx/clus_result.dat'
data = read_vmd_clus_result(path)
id_with_time = get_id_with_time(data)
plot_clustering_id_with_time(np.array(id_with_time)[:, 1], 0.5, path=os.path.dirname(path), point=data[0][0], ssize=1.25)

```

## FEP单点突变热图

读取数据：

```python
ddG = read_single()  
# not provided here. customize yourselves. it's just a dictionary of mutation: ddG.
# you must follow the format of E1A, E10A, etc.
ddG = {
	'E1A': -0.783225000000002,
 	'V2A': 0.379990000000001,
 	'T3A': -0.7186525,
 	'E4A': 2.6721,
    .....
}
```

然后进行热图绘制：

```python
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import copy
import seaborn as sns
# also requires the above plot_common

def get_matrix(ddG):
    columns = sorted(list(set([key[:-1] for key in ddG.keys()])), key=lambda x: int(x[1:]))
    rows = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z']
    df = pd.DataFrame(index=rows, columns=columns)
    df.iloc[:,:] = np.NAN
    # process into a matrix
    for key, value in ddG.items():
        col, row = key[:-1], key[-1]
        df.loc[row, col] = value
    return df
    
def heatmap_single(df):
    # Convert all entries to numeric values, replacing non-numeric entries with NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    # NOTE: you can adjust the color here. The effect may vary with your ddG data range.
    cmap = sns.diverging_palette(h_neg=0, h_pos=240, n=15, as_cmap=True)  
    # Mask for NaN values
    mask = df.isnull()
    # heatmap. no text means not done. the text color is adaptive to the background color
    fig, ax = plot_common('Residue Position', 'Mutant', size=(8, 10), xpad=6, ypad=0, ticks_size=14, tight=False)
    sns.heatmap(df, cmap=cmap, center=0, annot=True, mask=mask, fmt='.2f',
                cbar_kws={'label': '\u0394\u0394G (kcal/mol)', 'format': '%.2f'},
                linewidths=0.5, linecolor='grey',
                annot_kws={'fontfamily': 'Arial'})
    # Update colorbar font size
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('\u0394\u0394G (kcal/mol)', fontsize=18, weight='demibold', family='Arial')
    plt.show()


df = get_matrix(ddG)
heatmap_single(copy.deepcopy(df))
```