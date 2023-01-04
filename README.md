# low-genus-actions

This repository contains sage code to help determine if a group action on a surface is arithmetic.

Given a finite cover of a surface with deck group G:
* The module `surfaces` can construct the homology of the cover as a G-module.
* The module `lifted_twists` can check which Dehn twists lift and compute the action of (powers of) Dehn twists on the homology of the cover.
* The module `index` can determine whether a subgroup of SL(2,Z) has finite index.

Files starting with `2-` or `3-` are computations for specific cases.  They all show that the cover in question is arithmetic (running any of these files will print `True`).

Note that you should load the GAP package `FGA`.