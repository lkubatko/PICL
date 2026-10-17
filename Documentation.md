<h1>PICL: Phylogenetic Inference with Composite Likelihood</h1>

## Introduction

PICL stands for Phylogenetic Inference with Composite Likelihood. The PICL package includes a collection of methods for analyzing sequence data using composite likelihood.  It includes methods for estimating gene trees as well as those for estimating species trees under the multispecies coalescent (MSC). For analyses under the MSC, both multilocus and SNP data can be analyzed, and rate variation among sites can be incorporated. Multiple individuals per species can be included. Details about the method can be found in the following publication:

Kubatko, L., S. Kong, E. Webb, and Z. Chen. 2025. The promise of composite likelihood for species-level phylogenomic inference, <i>Evolutionary Journal of the Linnean Society</i> [<a href="https://doi.org/10.1093/evolinnean/kzaf008">web link</a>]

<br><br>

## Installation
PICL is written in C. After cloning the repository or downloading the code, PICL will need to be compiled.   The MacOS version of PICL uses routines from the <a href="https://github.com/boostorg/boost">Boost C++ libraries</a>. These will need to installed prior to compilation (see below for instructions).

Open a terminal window, navigate to the directory where the source code is stored, and issue the command: `make`. The program can then be called with `./picl`.

It is recommended that you run the two examples below before trying your own data. The first example will take about 30 minutes to run (reducing the number of optimization steps will reduce the time). The second should run in under 5 minutes. 

<br><br>

## Input file format
PICL requires two input files and will accept an optional third file. The first required file is the aligned sequence data, which should be formatted as a PHYLIP file. The loci/SNPs should be concatenated into one long alignment with no spaces (but please note that this is NOT a concatenation method). This file must be called `data.phy`.

The second required file is a settings file that provides instructions for what analyses the program should run, as well as information about assigning lineages to species.  The format is shown in the below. This file must be called `settings`. 



<details>
  <summary>Show settings file</summary>
Model: 1 <br>
Gaps: 1 <br>
Boot_type: 1<br>
Nboot: 100 <br> 
Theta: 0.004 <br>
Alpha: 10.0 <br>
Lambda: 500.0 <br>
Rate_param: 1.0<br>
Random_tree: 0<br>
Opt_bl: 1 <br>
User_bl: 1<br>
Num_opt: 100000 <br>
Seed1: 12345<br>
Seed2: 67890 <br>
Num_cat: 4 <br>
Tree_search: 0<br>
Num_iter: 10000 <br>
Multi_iter: 3<br>
Prob_bound: 0.05 <br>
Test_incr: 250 <br>
Opt_slope: -0.01 <br>
Beta: 0.005 <br>
Verbose: 1 <br>
Nspecies <br>
list of species names <br>
species_name lineage_name_1 <br>
species_name lineage_name_2 <br>
.
.
.
.
</details>

In the settings file, each keyword (i.e., the words followed by a colon) should appear exactly as written above, <b>in the same order,</b> and should be followed by a colon and then a space. Following the space, a number with no commas should appear. Below are the options for each keyword.

<ul>
<li> Model:
<ul>
<li> 1 = multilocus or CIS data
<li> 2 = multilocus or CIS data with discrete gamma-distributed rate variation 
<li> 3 = SNP data
<li> 4 = JC69 model for gene trees (this computes composite likelihood under the JC69 model without the multispecies coalescent)
<li>5 = multilocus or CIS data with variable population sizes
</li>
</ul>
<br>
<li> Gap:
<ul>
<li> 0 = sites containing a gap in at least one lineage are ignored
<li> 1 = all sites are included; gapped sites are currently ignored at the quartet level
</ul>
<br>
<li>Boot_type: if 0, no bootstrapping is done; if 1, tree bootstrapping is done; if 2, branch length bootstrapping on a fixed tree is done. In both cases, results are written to file boots.dat
<br>
  
<li> Nboot: number of bootstrap replicates for either bootstrap option
<br>


<li> Theta: a value for the effective population size parameter, $4N_e\mu$; this parameter will eventually be estimated, but a fixed value specified by the user is currently required
<br>

<li>Alpha: a value for the $\alpha$ parameter in the variable population size model
<br>

<li> Lambda: a value for the parameter $\lambda$ in the variable population size model
<br>

<li> Rate_param: value for the rate parameter in the discrete gamma model of rate variation; this will be optimized if branch lengths are optimized or a tree search is performed
<br>

<li> Random_tree: if 0, the tree will be read from treefile.tre; if 1, a random tree will be generated as a starting point for the search
<br>

<li> Opt_bl: 
<ul>
<li> 0 = no branch length optimization is done; the composite likelihood will be evaluated for the user-specified (if Random_tree = 0) or randomly generated (if Random_tree = 1) tree
<li> 1 = branch lengths on the input tree will be optimized using an uphill search; the number of iterations is given by Num_opt
<li> 2 = branch lengths on the input tree will be optimized using a simulated annealing search; the number of iterations is given by Num_opt
<li> 3 = branch lengths on the input tree will be optimized numerically using derivatives; this option is not yet implemented
</ul><br>

<li>User_bl: if 1, user branch lengths are used as the starting point for inference; if 0, branch lengths are ``evenly'' spaced throughout the tree. This option is ignored if Opt_bl=0.
<br>

<li> Num_opt: if Opt_bl = 1 or 2, this specifies the number of optimization steps
<br>

<li> Seed1: and Seed2: random number seeds (allows repeatability of analyses)
<br>

<li> Num_cat: if Model = 2, this specifies the number of categories for the discrete gamma model; it is ignored if Model $\neq$ 2
<br>

<li> Tree_search: specifies the type of tree search
<ul>
<li> 0 = none (user specified or randomly generated tree will be evaluated) 
<li> 1 = uphill search (NNI proposals)
<li> 2 = simulated annealing search (NNI proposals)
</ul>
<br>
<li> Max_iter: maximum number of iterations for the tree search (both methods)
<br>
<li> Multi_iter: number of iterations for tree scaling (both methods)
<br>

<li> Prob_bound: probability that a node has NOT been selected for rearrangement; this is used to specify a stopping rule in the annealing tree search and should be set to some small number (recommended is 0.05); this value is ignored unless Tree_search = 2
<br>

<li> Test_incr: specifies the frequency with which the cooling schedule  is updated in the adaptive annealing procedure; this value is ignored unless Tree_search = 2
<br>

<li> Opt_slope: specifies the target for adaptive annealing; this value should be set to -0.01 and should generally not be modified by users; this value is ignored unless Tree_search = 2
<br>

<li> Beta: parameter to control the cooling rate in the simulated annealing search (ignored if Tree_search = 0)
<br>

<li> Verbose: print extra checks and information? 0 = no, 1 = yes (option 1 used for de-bugging; users will want this set to 0)
</ul>

Following this, the number of species is listed. The next line contains a list of the species names. The remaining lines are used to assign lineages to species. The first word in each line must be one of the species names from the list provided. The second word in each line should be the lineage name <b>exactly as it appears in the PHYLIP input file (`data.phy`).</b> The lineages can appear in any order.

As the settings file is read, the information will be printed to the console. Users should check that all model options are being set as expected.  You should also check that the lineage-to-species assignments are being read correctly. 

The third, optional input file is a treefile formatted in Newick format that must be called `treefile.tre`. This should be provided if a user would like to estimate speciation times for a fixed species-level phylogeny. It can optionally be provided as a starting point for a tree search. Note that in version 1.0, branch lengths MUST be provided in the input tree (even if branch lengths will be optimized).

All files must be included in the directory from which PICL is called and must be named as above (i.e., `data.phy`, `settings`, and `treefile.tre`).

<br><br>


## Output files

Much of the information of interest to users is written directly to the console, so users should read through the output. In addition, optimized trees will be written to a file called `outtree.tre`.  If the analysis involves optimization along a fixed tree, this file will contain the tree with optimized branch lengths as a Newick string in both coalescent units and mutation units. If the analysis includes a tree search, then this file will contain both the current tree (i.e., the tree at the last iteration of the search) as well as the best tree encountered during the search in both coalescent units and mutation units, all in Newick format. 

If a bootstrapping analysis was carried out in order to estimate confidence intervals for speciation times on a fixed tree, the results of the bootstrap replicate analyses will be written to a file called `boots.dat`. The estimated speciation times for each bootstrap replicate are written in same order, with the times for one replicate per line (so the number of lines in the file will be equal to the number of bootstrap replicates). In order to determine the order in which the times are printed, the user must currently do a little work (improving the interface is on the list of future work) -- this is most easily done by comparing the printed times with the Newick tree that appears in `outtree.tre`.  Confidence intervals are obtained by importing the data into other software (e.g., R, Excel) and finding the percentiles to achieve the desired coverage. For example, for a 95\% confidence interval, one would find the $2.5^{th}$ percentile and the $97.5^{th}$ percentile of each column to obtain the lower and upper bounds, respectively, of the confidence interval.


<br><br>

## Example
One of the example files included in the examples directory corresponds to a simulated data set with 100,000 bp (100 bp/gene for 1,000 genes) for 10 species with 2 lineages per species (i.e., there are a total of 20 sequences in the aligned data set). The species are labeled $A$ - $J$, and the lineages are denoted $A.1$, $A.2$, $B.1$, $B.2$, etc. The `data.phy` file is formatted as follows, where the two numbers at the top of the file are the number of lineages (which is equal to the number of lines that follow) and the number of sites in the complete alignment). There should be a space separating these numbers; the data then start on the next line. Lineage names are no more than 10 characters long. The sequence data begin in column 11.

<details>
  <summary>Show data.phy</summary>
  <pre>
  20 100000
A.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
A.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
B.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
B.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
C.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
C.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
D.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
D.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
E.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
E.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
F.1       AGGTGCAGACGAAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
F.2       AGGTGCAGACGAAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
G.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
G.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
H.1       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
H.2       AGGTGCAGACGGAGCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
I.1       AGGTGCAGACGGATCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
I.2       AGGTGCAGACGGATCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
J.1       AGGTGCAGACGGATCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
J.2       AGGTGCAGACGGATCCTGCCCTAGCTAGAAGTAAAGGCACCGCTT...
</pre>
</details>

Several example analyses for these data are described below.

### Analysis 1: Simulated annealing search for the species tree
Model 1 is used to carry out a simulated annealing search for the tree with 30,000 iterations with $\beta = 0.005$ starting from a random tree.  The settings file for this analysis is shown below.

<details>
  <summary>Show settings file</summary>
Model: 1<br>
Gaps: 1<br>
Boot_type: 0<br>
Nboot: 0<br>
Theta: 0.004<br>
Alpha: 10.0<br>
Lambda: 500.0<br>
Rate_param: 1.0<br>
Random_tree: 1<br>
Opt_bl: 1<br>
User_bl: 0<br>
Num_opt: 100000<br>
Seed1: 12345<br>
Seed2: 67890<br>
Num_cat: 4<br>
Tree_search: 2<br>
Num_iter: 30000<br>
Multi_iter: 3<br>
Prob_bound: 0.05 <br>
Test_incr: 250 <br>
Opt_slope: -0.01  <br>
Beta: 0.005<br>
Verbose: 1<br>
10<br>
A B C D E F G H I J<br>
A A.1<br>
A A.2<br>
B B.1<br>
B B.2<br>
C C.1<br>
C C.2<br>
D D.1<br>
D D.2<br>
E E.1<br>
E E.2<br>
F F.1<br>
F F.2<br>
G G.1<br>
G G.2<br>
H H.1<br>
H H.2<br>
I I.1<br>
I I.2<br>
J J.1<br>
J J.2<br>
</details>

Once the analysis is complete, the results are written to the file `outtree.tre`:

<details>
<summary>Show outtree.tre</summary>
Current tree in mutation units: 
((A:0.005826,((D:0.001879,E:0.001879):0.001863,(B:0.001832,C:0.001832):0.001910):0.002084):0.001990,(F:0.005627,((I:0.001926,J:0.001926):0.002050,(G:0.001985,H:0.001985):0.001991):0.001651):0.002189);<br>

Current tree in coalescent units: 
((A:2.912863,((D:0.939681,E:0.939681):0.931364,(B:0.915871,C:0.915871):0.955173):1.041819):0.995043,(F:2.813526,((I:0.962947,J:0.962947):1.025063,(G:0.992345,H:0.992345):0.995665):0.825515):1.094381);<br>

Best tree found by the algorithm in mutation units: 
((A:0.005826,((D:0.001879,E:0.001879):0.001863,(B:0.001832,C:0.001832):0.001910):0.002084):0.001990,(F:0.005627,((I:0.001926,J:0.001926):0.002050,(G:0.001985,H:0.001985):0.001991):0.001651):0.002189);<br>

Best tree found by the algorithm in coalescent units: 
((A:2.912863,((D:0.939681,E:0.939681):0.931364,(B:0.915871,C:0.915871):0.955173):1.041819):0.995043,(F:2.813526,((I:0.962947,J:0.962947):1.025063,(G:0.992345,H:0.992345):0.995665):0.825515):1.094381);<br>
</details>

<br>
It is recommended that the search be run several times from multiple random starting trees. The tree with the highest composite likelihood across the replicates is the composite likelihood estimate. This can be done by running PICL with different random number seeds.
<br>

### Analysis 2:  Confidence intervals for speciation times along the user-provided species treee
In this analysis, we consider estimating speciation times by maximzing the composite likelihood, and we construct confidence intervals for the speciation times using the bootstrap. The tree for which the analysis should be carried out is provided in the file treefile.tre:

<details>
<summary>Show treefile.tre</summary>
((A:3.0,((B:1.0,C:1.0):1.0,(D:1.0,E:1.0):1.0):1.0):1.0,(F:3.0,((G:1.0,H:1.0):1.0,(I:1.0,J:1.0):1.0):1.0):1.0);
</details>

The settings file that corresponds to this analysis with 100 bootstrap replicates is shown below.

<details>
<summary>Show settings file</summary>
Model: 1 <br>
Gaps: 1 <br>
Boot_type: 2<br>
Nboot: 100<br>
Theta: 0.004 <br>
Alpha: 10.0<br>
Lambda: 500.0<br>
Rate_param: 1.0 <br>
Random_tree: 0 <br>
Opt_bl: 1 <br>
User_bl: 0 <br>
Num_opt: 1000000 <br>
Seed1: 12345 <br>
Seed2: 67890 <br>
Num_cat: 4 <br>
Tree_search: 0 <br>
Max_iter: 30000 <br>
Multi_iter: 3<br>
Prob_bound: 0.05 <br>
Test_incr: 250<br>
Opt_slope: -0.01 <br>
Beta: 0.005 <br>
Verbose: 1 <br>
10 <br>
A B C D E F G H I J <br>
A A.1 <br>
A A.2 <br>
B B.1 <br>
B B.2 <br>
C C.1 <br>
C C.2 <br>
D D.1 <br>
D D.2 <br>
E E.1 <br>
E E.2 <br>
F F.1 <br>
F F.2 <br>
G G.1 <br>
G G.2 <br>
H H.1 <br>
H H.2 <br>
I I.1 <br>
I I.2 <br>
J J.1 <br>
J J.2 <br>
</details>

As the analysis runs, the number of the current bootstrap replicate is written to the console. Upon completion, the optimized tree is written to `outtree.tre` and the speciation time estimates from the bootstrap replicates are written to the file `boots.dat`. Upon completion, that file should look like the following:

<details>
<summary>Show boots.dat</summary>
0.019415355440 0.011003931292 0.003590537588 0.013524495846 0.019415351969 0.011003981140 0.011010186935 0.013524499070 0.014524638664 <br>
0.017018017911 0.011534325721 0.003429927693 0.013519423921 0.017017926619 0.011534394471 0.012000000849 0.013519471464 0.014519483694 <br>
0.015253917685 0.011418249849 0.003305348187 0.013479000847 0.015253876771 0.011418560203 0.011450611066 0.013479239464 0.014479245549 <br>
0.016726396010 0.011999789116 0.003302720318 0.013149364175 0.016726306282 0.011999812719 0.012000000849 0.013149375931 0.014149389410 <br>
0.017899360190 0.011851918203 0.003437334294 0.013710358162 0.017899193012 0.011851965534 0.011876357250 0.013710425990 0.014710445530 <br>
0.016778388737 0.011950782721 0.003482994597 0.013834426784 0.016778317732 0.011950800529 0.012000000849 0.013834642007 0.014834677729 <br>
.
.
.
.
(100 lines total)
</details>

This file contains 9 columns 9 and 100 lines (they columns may have wrapped on the display). The 9 columns correspond to the 9 speciation times in the input tree. For example, the first column is the time of the root node, the second column is the time at which species $A$ speciates from species $B$, $C$, $D$, and $E$, etc. To figure out which times in the output correspond to which nodes in the tree, compare this file to the file `outtree.tre` -- the bootstrap estimates should be close to the estimates for the original data. The bootstrap data will need to be read into another software package (e.g., R, Excel) to be summarized.

<br><br>

## Algorithmic details 

<b> Random starting tree: </b>
If random_tree = 1, the method starts from a random tree generated from a Yule model. The starting branch lengths are "evenly" spaced throughout the tree, with the time of the root node set to be equal to the number of species divided by 2 (in coalescent units; this is equal to theta times the number of species divided by 2 in mutation units).<br>

<b> Theta:</b> The parameter $\theta$ is the effective population size parameter, which in PICL is $4N_e \mu$, where $N_e$ is the effective population size and $\mu$ is the mutation rate. Note that PAUP* uses $\theta/2$.<br>

<b>Rate variation with the discrete gamma model:</b> The current optimization alternates between branch length optimization and optimization of the parameter in the rate variation model. Specifically, Num_opt iterations are spent optimizing the branch lengths with the starting value for Rate_param fixed; then, Num_opt iterations are used to optimze the rate parameter with the branch lengths fixed; finally, Num_opt iterations are used for another round of branch length optimization with the rate parameter fixed. As a consequence, the number of total branch length iterations here is 2*Num_opt. A better procedure would be to integrate the optimization of all parameters (planned future work).

The maximum value for the rate parameter is set to 1000.0. This indicates that a model with no rate variation is likely appropriate. In this case, users may want to re-run their analysis with model=1 (instead of model=2), as the composite likelihood score may improve a bit without modeling rate variation (due to the upper limit on the rate parameter in model=2).

<b>Simulated annealing search:</b> The version of the simulated annealing search that is currently implemented is as follows:
<ul>
<li> A fixed number of iterations (1000) is used for the burnin period to set parameters in the cooling schedule. After burnin, an adaptive annealing schedule is used to tune the algorithm. The adaptive step is carried out every 250 iterations by default.
<li> The algorithm terminates when no new moves have been accepted after a specified number of proposals, and the chance of proposing a new accepted move is very low OR when the maximum user-specified number of iterations is reached. If the algorithm returns a number of iterations that is a multiple of 1,000, you should raise the maximum and re-run the algorithm.
<li> Following an NNI move, a new time for the affected node is generated uniformly over the allowable interval.  A better strategy is to use derivative information to propose a fairly good value for the new time.
<li> In the rate variation model, a new time for the node affected by the NNI is proposed jointly with a new value of the rate parameter. 
</ul>

<b>Bootstrapping:</b> The bootstrap procedure can be used to generate confidence intervals on the speciation times or to estimate support values for nodes. For each of the user-specific number of replicates, a bootstrap data set will be generated by sampling sites from the original alignment to generate a new alignment. 

<p> If boot_type = 2, speciation times will be estimated from the new alignment, and will be written to a file called boots.dat. Upon completion of the user-specified number of bootstrap replicates, the columnus in the boots.dat file can be summarized to obtain the endpoints of the confidence interval. Users will need to ``match'' the numerical values printed to the tree in the file outtree.tre to understand the ordering of the speciation times. </p>

<p>If boot_type=1, the simulated annealing algorithm will be used to estimate a tree from each bootstrap sample. Results for each replicate will be written to the file boots.dat. These trees will need to be summarized in another software (e.g., PAUP*) to obtain the support values on the nodes.</p>

<p>Note that the multilocus bootstrap is not implemented. </p>



<br><br>

## Future work and software comparison

The composite likelihood calculation implemented here matches that in PhyNEST (Kong et al., Syst. Biol.  74(1): 53–69, 2025) when the number of reticulations is set to 0 and the data do not include gaps. It does not match that implemented in qAge in PAUP* because the PAUP* implementation further reduced the site pattern counts down to 9 or 11 categories to increase computational efficiency on fixed trees (see paper referenced above for details).

Several extensions of this method are under active development in the Kubatko Research Group, with planned incorporation in this software. In particular, methods for testing hypotheses using composite likelihood ratio tests are being developed and a Bayesian species tree inference method using composite likelihood is in the works. In addition, the following more minor additions to the software are planned:

<ul>
<li> Remove requirement for the tree in the file treefile.tre to have branch lengths</li>
<li> Bootstrap support values on nodes [done, v. 1.0.5+, 7/22/26]
<li> Additional options for handling gaps and ambiguity codes (currently ignored at the quartet level) [ambiguity codes are supported in v 1.0 and above]
<li> Derivative-based branch length estimation
<li> Efficient calculation of CL after NNI 
<li> Parallel computation of CL
<li> Estimation of the effective population size parameter $\theta$ (models 1, 2, and 3)
<li> Posterior site pattern probabilities (Richards and Kubatko, 2021)
<li> Site pattern probabilities under a relaxed clock (Richards and Kubatko, 2022)
<li> Site pattern probabilities for gene tree estimation [done, 7/20/25]
</ul>


## Installation information for Boost libraries (MacOS)

The <a href="https://github.com/boostorg/boost">Boost libraries</a> can be installed with either of these two methods:
<ul>
<li> <b>Using homebrew:</b> `brew install boost`
<li> <b>Manual install:</b> Go to <a href="https://www.boost.org/users/download/">https://www.boost.org/users/download/</a> and download the appropriate .gzip file. Expand the file (e.g., `gunzip boost_1_91_0.tar.gz`) and move the tar file into the directory from which you'll run PICL. Then issue the command `tar xzf boost_1_91_0.tar`.
</ul>
