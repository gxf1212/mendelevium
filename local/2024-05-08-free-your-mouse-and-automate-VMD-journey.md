---
layout: post
title: "Free your mouse and automate your VMD journey!"
categories: "Visualization"
tags: MD VMD
author: Xufan
---

# Free your mouse and automate your VMD journey!

VMD provides a lot of useful features. Though the GUI is ugly, VMD is irreplacable, especially when dealing with trajectories. 

You must have heard that "we can find an equivalent command for everything you can do through clicking". This is usually true for VMD, Pymol and UCSF Chimera. You always need to check a few similar trajectories that you hope to process similarly, e.g. replicas, or control groups. You want to save all the configuration code and reuse them.

VMD scripting is powerful, too. But plugins may have packaged a few functions

## Basic workflow

### Global setting

```tcl
display projection Orthographic
color Display {Background} white
color Labels {Atoms} black
color Labels {Bonds} mauve
label textthickness 2.5
```

We don't use "Perspective", 失真

### Load data

```tcl
mol new npt.gro type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all  # add molecule
mol addfile md_pbcmol.xtc type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all # load trajectory
animate delete beg 0 end 0  # remove the first frame, if you use .gro/.pdb file
set numframes [molinfo top get numframes]  # get number of frames
```



> [!TIP]
>
> We assume you have done PBC processing. VMD provides PBC Tools, which don't always work. Follow us for later posts about this topic!

### Add molecule representation



```tcl
mol representation CPK 0.200000
mol color Type
mol selection {(protein or water) and (same residue as within 5 of residue LIG)}
mol material Opaque
mol addrep top
# update selection/color every frame, 1 means yes
mol selupdate 3 top 1
mol colupdate 3 top 0
```



> [!NOTE]
>
> Replacing `LIG` with `$ligname` (so that you can dynamically assign) does not work...so try to use the same residue name for different ligands.



### Save visualization state



## What does a plugin look like





how to envoke a plugin

RMSD TT





## End

In another article, we will introduce examples of calling plugins from cmd





# Three pieces of code that help you automate post-analysis of protein-ligand system in VMD



## RMSD Trajectory Tool





## Clustering analysis

### Algorithm





### Plotting the result



```python

```





## Movie Maker




