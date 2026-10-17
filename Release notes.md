<h1>Release notes</h1>

<h2>Version 1.2 - October 17, 2026</h2>
<ul>
<li> Corrected a bug in printing branch lengths for bootstrap replicates; branch lengths are printed in mutation units
<li> Improved Makefile to optimize compilation and to fix reproducibility errors due to floating point rounding
<li> Improved Makefile to utilize boost libraries on unix/linux HPC clusters
<li> Implemented caching of quartet site pattern counts in composite likelihood computation to speed computation
<li> Improved algorithm for proposing a new branch length during simulated annealing
<li> Added checks for keyword ordering in the settings file
<li> Added checks for matches between lineage names in settings file and lineage names in data.phy; removed restrictions on ordering in these files
<li> Previously, PICL read in theta = 2Nu (even though the documentation stated 4Nu); this is now corrected: PICL reads in theta=4Nu from the settings file
<li> The settings file has been updated to include an improved version of model 5 
</ul> 

<ul><li>Known bugs:<br>
<ul>
<li> Bootstrapping occasionally gives an error due to the interface between PICL and the random number generation routines.
</ul>
</ul>

<h2>Version 1.0.5 - July 22, 2026</h2>
<ul>
  <li>Added bootstrap across trees:  new option in settings file Boot_type is set to 1 for tree bootstrap, and to 2 for branch length bootstrap on a fixed tree</li>
</ul>

