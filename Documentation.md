<h1>PICL: Phylogenetic Inference with Composite Likelihood</h1>

## Introduction

PICL stands for Phylogenetic Inference with Composite Likelihood. The PICL package includes a collection of methods for analyzing sequence data using composite likelihood.  It includes methods for estimating gene trees as well as those for estimating species trees under the multispecies coalescent (MSC). For analyses under the MSC, both multilocus and SNP data can be analyzed, and rate variation among sites can be incorporated. Multiple individuals per species can be included. Details about the method can be found in the following publication:

Kubatko, L., S. Kong, E. Webb, and Z. Chen. 2025. The promise of composite likelihood for species-level phylogenomic inference, <i>Evolutionary Journal of the Linnean Society</i> [<a href="https://doi.org/10.1093/evolinnean/kzaf008">web link</a>]

<br><br>

## Installation
PICL is written in C. After cloning the repository or downloading the code, PICL will need to be compiled. 

On a Mac, open a terminal window, navigate to the directory where the source code is stored, and issue the command: `gcc main.c -lm -o picl`. The program can then be called with `./picl`.

It is recommended that you run the two examples below before trying your own data. The first example will take about 30 minutes to run (reducing the number of optimization steps will reduce the time). The second should run in under 5 minutes. 

<br><br>

## Input file format
PICL requires two input files and will accept an optional third file. The first required file is the aligned sequence data, which should be formatted as a PHYLIP file. The loci/SNPs should be concatenated into one long alignment with no spaces (but please note that this is NOT a concatenation method). This file must be called `data.phy`.

The second required file is a settings file that provides instructions for what analyses the program should run, as well as information about assigning lineages to species.  The format is shown in the below. This file must be called `settings`. 



<details>
  <summary>Show settings file</summary>
Model: 1 <br>
Gaps: 1 <br>
Bootstrap: 0 <br> 
Theta: 0.002 <br>
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
</ul>
<br>
<li> Gap:
<ul>
<li> 0 = sites containing a gap in at least one lineage are ignored
<li> 1 = all sites are included; gapped sites are currently ignored at the quartet level
</ul>
<br>
<li> Bootstrap: if 0, no bootstrapping is done; if >0, this specifies the number of bootstrap replicates used to obtain confidence intervals on the speciation times of a fixed species tree. The Opt_bl option below must not be set to 0 to carry out bootstrapping.
<br>

<li> Theta: a value for the effective population size parameter; this parameter will eventually be estimated, but a fixed value specified by the user is currently required
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
<li> Num_iter: number of iterations for the tree search (both methods)
<br>

<li> Beta: parameter to control the cooling rate in the simulated annealing search (ignored if Tree_search = 0)
<br>

<li> Verbose: print extra checks and information? 0 = no, 1 = yes (option 1 used for de-bugging; users will want this set to 0)
</ul>

Following this, the number of species is listed. The next line contains a list of the species names. The remaining lines are used to assign lineages to species. The first word in each line must be one of the species names from the list provided. The second word in each line should be the lineage name <b>exactly as it appears in the PHYLIP input file (`data.phy`).</b> The lineages can appear in any order.

As the settings file is read, the information will be printed to the console. Users should check that all model options are being set as expected.  You should also check that the lineage-to-species assignments are being read correctly. 

The third, optional input file is a treefile formatted in Newick format that must be called `treefile.tre`. This should be provided if a user would like to estimate speciation times for a fixed species-level phylogeny. It can optionally be provided as a starting point for a tree search.

All files must be included in the directory from which PICL is called and must be named as above (i.e., `data.phy`, `settings`, and `treefile.tre`).

<br><br>


## Output files

Much of the information of interest to users is written directly to the console, so users should read through the output. In addition, optimized trees will be written to a file called `outtree.tre`.  If the analysis involves optimization along a fixed tree, this file will contain the tree with optimized branch lengths as a Newick string in both coalescent units and mutation units. If the analysis includes a tree search, then this file will contain both the current tree (i.e., the tree at the last iteration of the search) as well as the best tree encountered during the search in both coalescent units and mutation units, all in Newick format. 

If a bootstrapping analysis was carried out in order to estimate confidence intervals for speciation times on a fixed tree, the results of the bootstrap replicate analyses will be written to a file called `boot.dat`. The estimated speciation times for each bootstrap replicate are written in same order, with the times for one replicate per line (so the number of lines in the file will be equal to the number of bootstrap replicates). In order to determine the order in which the times are printed, the user must currently do a little work (improving the interface is on the list of future work) -- this is most easily done by comparing the printed times with the Newick tree that appears in `outtree.tre`.  Confidence intervals are obtained by importing the data into other software (e.g., R, Excel) and finding the percentiles to achieve the desired coverage. For example, for a 95\% confidence interval, one would find the $2.5^{th}$ percentile and the $97.5^{th}$ percentile of each column to obtain the lower and upper bounds, respectively, of the confidence interval.


<br><br>

## Example
The example files included in the examples directory correspond to a simulated data set with 100,000 bp (100 bp/gene for 1,000 genes) for 10 species with 2 lineages per species (i.e., there are a total of 20 sequences in the aligned data set). The species are labeled $A$ - $J$, and the lineages are denoted $A.1$, $A.2$, $B.1$, $B.2$, etc. The `data.phy` file is formatted as follows, where the two numbers at the top of the file are the number of lineages (which is equal to the number of lines that follow) and the number of sites in the complete alignment). There should be a space separating these numbers; the data then start on the next line. Lineage names are no more than 10 characters long. The sequence data begin in column 11.

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
Bootstrap: 0<br>
Theta: 0.002<br>
Rate_param: 1.0<br>
Random_tree: 1<br>
Opt_bl: 1<br>
User_bl: 1<br>
Num_opt: 100000<br>
Seed1: 12345<br>
Seed2: 67890<br>
Num_cat: 4<br>
Tree_search: 2<br>
Num_iter: 30000<br>
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
Current tree in mutation units: <br>
((A:0.005826,((B:0.001832,C:0.001832):0.001910,(D:0.001879,E:0.001879):0.001863):0.002084):0.001990,(F:0.005627,((G:0.001985,H:0.001985):0.001991,(I:0.001926,J:0.001926):0.002050):0.001651):0.002189);<br>
 <br>
Current tree in coalescent units: <br>
((A:2.912867,((B:0.915884,C:0.915884):0.955149,(D:0.939662,E:0.939662):0.931371):1.041834):0.995068,(F:2.813559,((G:0.992357,H:0.992357):0.995680,(I:0.962967,J:0.962967):1.025070):0.825522):1.094376); <br>
 <br>
Best tree found by the algorithm in mutation units:  <br>
((A:0.006188,((C:0.001352,B:0.001352):0.002746,(E:0.002185,D:0.002185):0.001913):0.002090):0.001456,(F:0.005149,((J:0.001066,I:0.001066):0.002587,(H:0.002019,G:0.002019):0.001634):0.001497):0.002494); <br>
 <br>
Best tree found by the algorithm in coalescent units:  <br>
((A:3.093902,((C:0.676158,B:0.676158):1.372808,(E:1.092525,D:1.092525):0.956441):1.044936):0.727879,(F:2.574742,((J:0.533016,I:0.533016):1.293316,(H:1.009541,G:1.009541):0.816791):0.748410):1.247039); <br>
</details>

<br>

### Analysis 2:  Confidence intervals for speciation times along the user-provided species treee
In this analysis, we consider estimating speciation times by maximzing the composite likelihood, and we construct confidence intervals for the speciation times using the bootstrap. The tree for which the analysis should be carried out is provided in the file treefile.tre:

<details>
<summary>Show treefile.tre</summary>
((A:4.0,((B:1.0,C:1.0):1.0,(D:1.0,E:1.0):1.0):1.0):1.0,(F:3.0,((G:1.0,H:1.0):1.0,(I:1.0,J:1.0):1.0):1.0):1.0);
</details>

The settings file that corresponds to this analysis with 100 bootstrap replicates is shown below.

<details>
<summary>Show settings file</summary>
Model: 1 <br>
Gaps: 1 <br>
Bootstrap: 100 <br>
Theta: 0.002 <br>
Rate_param: 1.0 <br>
Random_tree: 0 <br>
Opt_bl: 1 <br>
User_bl: 0 <br>
Num_opt: 1000000 <br>
Seed1: 12345 <br>
Seed2: 67890 <br>
Num_cat: 4 <br>
Tree_search: 0 <br>
Num_iter: 30000 <br>
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

As the analysis runs, the number of the current bootstrap replicate is written to the console. Upon completion, the optimized tree is written to `outtree.tre` and the speciation time estimates from the bootstrap replicates are written to the file `boot.dat`. Upon completion, that file should look like the following:

<details>
<summary>Show boot.dat</summary>
3.786275575362 2.804381515550 1.802258359673 0.806892917114 0.857299662938 2.790955951219 1.935663781232 0.935083170130 1.026693794973 <br>
3.976534692093 2.912370382952 1.906134787601 0.914466815066 0.936859157275 2.904599580271 2.053118571303 0.901589707544 1.032576239472 <br>
3.954904633952 3.015915360466 1.954042764148 0.949460451344 0.959741073899 2.894189894672 2.048634974469 1.162039017279 1.021140136548 <br>
3.969057565869 2.897533805650 1.901667393982 0.983553131189 0.928252523368 2.923025310521 2.013753600804 1.048998435202 0.939770528940 <br>
3.850675958129 2.868895441101 1.867414573296 0.855517748873 0.901075788504 2.777246709978 1.977418579396 0.977708280775 1.050581928432 <br>
3.874603336188 2.844417645252 1.818412711290 0.858541215538 0.828105209576 2.746019169706 1.991446075304 0.968493549217 1.004227258738 <br>
3.853837187110 2.777340087119 1.798691624615 0.886873902377 0.899984571724 2.866105901915 2.053879099898 0.967700687927 0.953500166143 <br>
3.933836477515 2.994085845869 1.895162004595 0.788799264008 1.026094394604 2.800655410989 2.021646178919 0.990970635605 0.998413233854 <br>
3.873853140418 2.858457177586 1.877808037618 0.897690275846 0.851489234741 2.845321219830 1.935706187312 0.939938218293 0.981964021395 <br>
3.947552189299 2.948826549072 1.860734756176 0.894474447332 0.965053036945 2.764652303285 1.957404678984 0.941276338507 0.954315092433 <br>
3.777940893859 2.850847482533 1.726640277141 0.914194393367 0.900477893534 2.759881119355 1.913459287407 0.926027611280 0.938402094386 <br>
3.831318631138 2.767511539937 1.762861295659 0.786982290692 0.894378585439 2.766984440469 1.908049119667 0.931050204043 0.874334311319 <br>
3.888243742647 2.913218374389 1.920349524830 0.955012308509 1.071823730611 2.744224886435 1.994352900152 0.974414052414 0.971082704540 <br>
3.905916291884 2.847767836403 1.876560034994 0.968354744578 0.928707302093 2.900358841208 2.016848865819 0.981445985513 0.908445159045 <br>
3.889252441324 2.885929171839 1.861894417453 0.867574197091 0.987221651421 2.832142404644 1.980005303413 1.039149529952 0.941672886994 <br>
3.818645593244 2.731056739323 1.782197644054 0.723256681049 0.919270715623 2.920263183515 2.139239114531 1.046991671035 1.005760196434 <br>
3.939164784335 2.918629295297 1.897821342376 0.970556156125 0.919656224308 2.760227511969 1.995366822710 1.015432991776 0.977081517129 <br>
3.949749501766 3.011063178587 1.897598250846 0.908144596216 0.929804724787 2.775693329705 1.929730574027 1.093901942040 0.980371295087 <br>
4.030247756398 2.913670210028 1.822157722992 0.952764458671 0.925002375352 2.891393262309 2.051403095568 1.025509578934 1.000851119028 <br>
3.807201848676 2.772934825956 1.768876140869 0.957957173218 0.859003227140 2.788236373874 1.928162719853 0.993781554652 0.994888669862 <br>
3.803806706916 2.790692580660 1.783075322235 0.822963212238 0.856451958259 2.815765484131 1.932115615776 1.045683789190 1.008707527685 <br>
3.808525837734 2.839942081577 1.824031541446 0.882837854616 0.903333130879 2.672621066716 1.866745384664 0.927628319228 0.864769430930 <br>
3.985192228570 2.923198538115 1.901130811322 0.928418825661 1.054115133396 2.816253330190 2.114601386583 1.013609409096 1.098141077267 <br>
3.896858476629 2.920013245050 1.884328819163 0.907646612673 1.005151307532 2.777537941195 1.952510471964 0.917100313775 0.958122158816 <br>
4.018059845174 2.996806998091 1.907223676931 0.894071917035 0.916257700422 2.878914854828 2.021834107643 1.000128031789 1.023210591736 <br>
3.795215096135 2.807908784274 1.803743134989 0.884664303153 0.862523521247 2.751476879013 1.918967213260 1.074820742338 0.938409361569 <br>
3.948318333166 2.944504682265 1.901076068715 0.888092432252 0.858821753072 2.798267164273 2.002126887741 0.971933458374 1.045773319620 <br>
3.995149386142 3.000482638568 1.920599619594 0.951488026279 0.950661184624 2.736098467043 2.057584142507 0.945392328437 0.994172563052 <br>
3.921564737843 2.957938288023 1.970120387557 0.977762904038 0.966932532125 2.810276460601 2.046349126224 1.026825182461 1.035209043856 <br>
3.934294491066 2.872468527457 1.897668509793 0.926489572226 0.982719035210 2.814890654883 2.057990784733 0.975855256526 1.013318326198 <br>
4.033746851085 2.965599609767 1.951049262025 1.002080905461 0.976214662391 2.841746870548 2.088183151348 1.027323928687 0.981688365437 <br>
3.985540994306 2.878336320682 1.824399262690 0.800815547539 0.843196643880 2.836539471253 1.946748426103 0.944977486886 0.973812594220 <br>
3.858706859532 2.980636627612 1.841619785624 0.948887014869 0.918737898646 2.765354411876 1.956748063174 0.989775491100 0.975878472234 <br>
4.004756822806 3.078531146051 1.931648756076 0.931926031014 0.939847105114 2.760007251736 2.010558011637 1.024848648559 1.032518917901 <br>
3.909899145943 2.953912128128 1.921618150698 0.955094977902 0.962934774880 2.725177847438 1.926100423682 0.917291077240 0.987267768555 <br>
3.951253218900 2.869230551799 1.767848347681 0.869012933286 0.839981230947 2.786157409582 2.030254556086 1.017381878634 0.992102754890 <br>
3.834313041178 2.873427607365 1.887162467983 0.905760969327 0.986325047788 2.814193831629 1.994935228774 0.984276708166 1.017815523980 <br>
3.947616643365 2.944828732624 1.892681365721 0.847614258524 0.896426973405 2.742164896092 2.008994000060 0.995164045602 0.947309543298 <br>
3.945489122657 2.936851223498 1.877964208966 1.013461673999 0.988444038780 2.931086637263 2.058069106566 1.039603502670 1.005775799728 <br>
4.018044111212 2.954780446771 1.853470737441 0.864779733722 0.988807665810 2.879256455223 2.027152025861 1.061484367086 1.046738216173 <br>
3.779521708844 2.804086136073 1.814696011544 0.905266585574 0.896158019397 2.715275870464 1.972784640331 0.976541100551 0.943339968718 <br>
3.971274105049 3.153092785406 1.966244851303 0.979341591011 1.002850991916 2.811848242716 1.983490089822 0.995098924645 0.945029106053 <br>
3.837282691827 2.944136933945 1.819762641463 0.868665480349 0.831586707570 2.690019654478 1.922082302096 0.996752675697 0.962867667001 <br>
3.869363235458 2.875149609784 1.916983379481 0.965284560877 0.914767364060 2.624823872048 1.837995317466 0.957112482247 0.889797222325 <br>
3.858619181527 2.769496040685 1.805171396740 0.872649419180 0.905381664304 2.842912029477 2.015085136782 1.007590161974 0.942310347191 <br>
3.975507675360 2.792084426739 1.835895474871 1.004211855704 1.006987774664 2.889165638115 2.033668343314 0.916987020928 0.923929697014 <br>
3.918223083375 2.874410275434 1.782732130243 0.861773273403 0.913391348293 2.722006122420 1.990047091857 0.959873096921 1.031290047367 <br>
3.811498967763 2.868679115764 1.922561067756 1.039233497093 0.918093323581 2.826608522962 1.982947566883 1.037805867264 0.981747092672 <br>
3.884446885012 2.898290548912 1.814698574568 0.796254215387 0.847578848418 2.817064041297 1.945009410651 0.978624715148 0.945133429642 <br>
3.876215826804 2.887347824861 1.913541054748 0.994748758498 0.928497984144 2.857343603446 1.931738118302 0.992832544346 0.965048182598 <br>
3.923485069722 2.973144391409 2.034412486626 1.005488908245 1.077776921236 2.851188679766 1.922266257784 0.984872702686 0.943553218010 <br>
3.936282369292 2.986620683536 1.885612331185 0.971258890392 0.918212943996 2.693471633109 1.966207567231 1.007359537012 0.938634453140 <br>
3.796487936033 2.904745590915 1.901567659040 0.915333826019 0.900368736753 2.733186524204 2.069066504086 1.041152403719 1.009015752533 <br>
3.917209828822 2.895341081616 1.760003428713 0.791878070509 0.928446137041 2.894133070882 2.117961069574 1.097771344493 1.025371241845 <br>
4.065724239767 2.972840163155 1.934527293025 0.981674461276 0.898469757499 2.911654943917 2.059366771050 1.029864684062 0.977771673300 <br>
3.997904965014 2.891944989528 1.877347956897 0.859145460750 1.007864991882 2.850566628870 1.977154163018 0.937055176905 0.916987019042 <br>
3.874586699491 2.954351297810 1.950984091330 0.975988710491 0.980846080336 2.742279543918 1.913072961511 0.920404531555 0.889859670965 <br>
4.005325133846 2.949239663438 1.890680007720 0.930667616600 0.989692597280 2.934621167218 2.045875895923 0.999903379769 1.021165355278 <br>
3.931828023525 2.871767651372 1.873346337540 0.851161095594 0.958313574442 2.959683162639 2.016297492457 1.027981445372 1.040596572649 <br>
3.913133422096 2.820329937687 1.881925858773 0.841471933476 0.865277333263 2.938636859902 2.011959691729 1.041902883485 0.988526557095 <br>
3.883631555450 2.918911127347 1.841488822741 0.902369648373 0.924596716724 2.891022885139 1.931531765846 0.858237637854 1.067064825406 <br>
3.924394132301 2.918263957064 1.960183215696 0.793629491641 0.899748088983 2.783717788335 1.921117787068 0.991697225351 0.965617951303 <br>
4.099589413102 3.041895840058 1.904848822939 0.916665219402 0.888296002547 2.894158668231 1.983017502702 1.030022906751 0.953579907762 <br>
4.086820556243 2.990635422036 1.889522966758 0.876558150835 0.931654937827 2.944680079319 1.995484639997 1.008556120840 0.971223117394 <br>
3.933660874654 2.745106577918 1.809256684898 0.800406903725 0.849854775998 3.104376212472 2.191916406731 1.021126001419 1.138052333569 <br>
3.907367754125 2.965789361509 1.907925951311 0.954659943189 0.897162312594 2.833080035952 1.972846002277 0.993648837296 0.943073746142 <br>
3.903432863522 2.978039103408 1.927181007851 1.029777923075 1.048162642341 2.833785031353 2.023465002735 1.079735641906 0.941335310829 <br>
3.945943947675 2.984181498837 1.921438007133 0.930825319900 1.054469344713 2.816046942776 1.971186392955 0.969373683751 1.001030174877 <br>
3.818710864437 2.805776915870 1.802503979458 0.847192822682 0.858160537319 2.700999225520 1.928498225059 0.982696339361 0.970140514269 <br>
3.946850616127 3.010264333777 2.056901454382 0.940613845442 0.987345945830 2.791213511972 1.916904616980 1.014488797791 0.909738455940 <br>
3.882424583439 2.873696747079 1.875328901143 0.937949038156 0.989516397958 2.841718392831 1.971857846627 0.995032212173 0.972875515839 <br>
3.811522026931 2.877611524666 1.821722153420 0.911141439750 0.973973418632 2.833345491247 1.916850958438 1.080285679103 0.815050817145 <br>
3.986452943005 2.908221233255 1.946772232562 0.949928259640 0.951791924555 2.906492574862 1.995931722414 0.995844385723 1.027652479906 <br>
3.877015168309 2.931455330782 1.884250980906 0.926514132590 0.880590618071 2.787531369186 1.921421975556 0.912563993456 1.016960158193 <br>
3.833067555092 2.853306693402 1.871496948076 0.919342859280 0.955475741111 2.799281187526 1.910706334349 0.936355151680 0.894651980358 <br>
3.914146578647 2.980585424527 1.869412445740 0.957110475995 0.904148798648 2.743448776932 1.924275041584 0.953480964329 0.970040286001 <br>
4.025901938707 2.952309550598 1.881142242517 0.915983897399 0.927976578546 2.828766744249 1.989617134408 0.921635371027 0.948510899666 <br>
3.979228226492 2.992259647240 1.935088636593 0.983974193850 1.071883318721 2.829081863914 1.937527713934 0.987257254360 0.921415147742 <br>
3.939564647499 2.933454083911 1.888771951983 0.940242540793 0.821425938152 2.999518596363 2.137269262276 0.933897576815 1.074751762293 <br>
4.067113789001 2.970581846953 1.884600611974 0.996032784896 0.938676101655 2.955264833864 2.030119958377 1.045563734141 0.986316455581 <br>
3.886135162601 2.889096683583 1.848660566548 0.825980206880 0.921633413462 2.816276425404 1.922463389160 1.095243008897 0.900182608514 <br>
3.957784549861 2.936493541488 1.819539981022 0.914363388395 0.940679125663 2.883856809160 2.098693422352 1.003347426191 1.046029962632 <br>
4.025939826900 3.011760204105 1.894213582818 0.956099927964 1.049888430404 2.881599202934 2.030560049513 1.047836344520 0.956383225077 <br>
3.961847850259 3.043898667913 1.894667831334 0.990448601082 1.040481115838 2.766956786247 1.999845908043 1.028857828531 0.896746612554 <br>
3.900760060882 2.827141463661 1.893337359171 0.978483163208 0.932324212757 2.850915243684 1.982896357014 0.972748537701 0.952585529479 <br>
4.076035348769 3.064217505073 1.903664134846 0.928866414646 0.980997682816 2.944374585253 2.024494282761 1.083482189163 1.037360948977 <br>
3.853062117473 2.893040092253 1.844360479673 0.878565859001 0.923381691596 2.840239706170 1.929725734005 0.910683613133 0.967751832058 <br>
3.840854253589 2.911735104734 1.734679281174 0.779758644482 0.821284172927 2.842115346140 1.998591741153 0.982635653737 1.010360136694 <br>
3.832495811728 2.970496750236 1.928169321301 0.907113194961 1.018047422437 2.810946677375 2.016752670853 1.090045172159 0.939641560312 <br>
3.859480914444 2.910132116218 1.872529094163 0.928774080691 1.010360482334 2.692848312010 1.891345275860 0.916884876612 0.901959416808 <br>
3.868600503925 2.984665244766 1.880777225831 0.857483100896 0.915886191379 2.767087711054 1.926844733069 1.031319218040 0.883411963770 <br>
3.850166869919 2.822736447763 1.836756763355 0.852832358811 1.001248883981 2.845908251999 1.931316231578 0.925963442627 0.970574721130 <br>
3.927589742424 2.869454735348 1.898494706529 0.944939701864 0.969295290019 2.799301884481 2.029333705108 1.080028772805 0.954948637831 <br>
4.028765688318 3.001831117685 1.895270280121 0.923952474254 1.035145523121 2.908831470549 1.981986041519 0.971336638410 0.921196890032 <br>
3.791795641648 2.958410363597 1.882022991275 0.936996837912 0.887934239647 2.650734900927 1.886587823018 0.958523225827 0.944934998009 <br>
3.844356485441 2.940268128247 1.858430601441 0.946670046679 0.976599570228 2.800917838788 2.069428884514 1.001801730647 1.090553925576 <br>
3.784481217630 2.934887484802 1.924346897825 0.971307355069 0.968814006181 2.747230814503 1.932851485829 0.952609901819 0.909690354166 <br>
3.854429354927 2.866044191654 1.784169431144 0.873865433299 0.876977378530 2.708820650240 1.910127795007 0.987685511042 0.971993628216 <br>
3.779462857146 2.928840024467 1.916813245762 0.947443832827 0.997932095963 2.647967731950 1.918345711809 0.961460699910 0.918069881490 <br>
3.759965819840 2.903053200788 1.869472655211 0.897924016642 0.924124285965 2.713976390169 1.979991395156 1.030592415149 0.880743556399 <br>
</details>

This file contains 9 columns 9 and 100 lines (they columns may have wrapped on the display). The 9 columsn correspond to the 9 speciation times in the input tree. For example, the first column is the time of the root node, the second column is the time at which species $A$ speciates from species $B$, $C$, $D$, and $E$, etc. To figure out which times in the output correspond to which nodes in the tree, compare this file to the file `outtree.tre` -- the bootstrap estimates should be close to the estimates for the original data. The bootstrap data will need to be read into another software package (e.g., R, Excel) to be summarized.

<br><br>

## Algorithmic details 

<b> Random starting tree: </b>
If random_tree = 1, the method starts from a random tree generated from a Yule model. The starting branch lengths are "evenly" spaced throughout the tree, with the time of the root node set to be equal to the number of species divided by 2 (in coalescent units; this is equal to theta times the number of species divided by 2 in mutation units).<br>

<b> Theta:</b> The parameter $\theta$ is the effective population size parameter, which in PICL is $4N_e \mu$, where $N_e$ is the effective population size and $\mu$ is the mutation rate. Note that PAUP* uses $\theta/2$.<br>

<b>Rate variation with the discrete gamma model:</b> The current optimization alternates between branch length optimization and optimization of the parameter in the rate variation model. Specifically, Num_opt iterations are spent optimizing the branch lengths with the starting value for Rate_param fixed; then, Num_opt iterations are used to optimze the rate parameter with the branch lengths fixed; finally, Num_opt iterations are used for another round of branch length optimization with the rate parameter fixed. As a consequence, the number of total branch length iterations here is 2*Num_opt. A better procedure would be to integrate the optimization of all parameters (planned future work).

The maximum value for the rate parameter is set to 1000.0. This indicates that a model with no rate variation is likely appropriate. In this case, users may want to re-run their analysis with model=1 (instead of model=2), as the composite likelihood score may improve a bit without modeling rate variation (due to the upper limit on the rate parameter in model=2).

<b>Simulated annealing search:</b> The version of the simulated annealing search that is currently implemented is naive in several ways:
<ul>
<li> A fixed number of iterations (100) is used for the burnin period to set parameters in the cooling schedule. The cooling schedule is not currently optimized, and problem-specific information is not used for tuning.
<li> Following an NNI move, a new time for the affected node is generated uniformly over the allowable interval.  A better strategy is to use derivative information to propose a fairly good value for the new time.
<li> No stopping rule is currently implemented, so the algorithm will use the maximum number of iterations specified by the user, even in the case that most proposed moves are being rejected.
<li> In the rate variation model, a new time for the node affected by the NNI is proposed jointly with a new value of the rate parameter $\alpha$. The 
</ul>

<b>Bootstrapping:</b> The bootstrap procedure can be used to generate confidence intervals on the speciation times. For each of the user-specific number of replicates, a bootstrap data set will be generated by sampling sites from the original alignment to generate a new alignment. Speciation times will be estimated from the new alignment, and will be written to a file called boot.dat. Upon completion of the user-specified number of bootstrap replicates, the columnus in the boot.dat file can be summarized to obtain the endpoints of the confidence interval. Users will need to ``match'' the numerical values printed to the tree in the file outtree.tre to understand the ordering of the speciation times. 

Note that the multilocus bootstrap is not implemented. Note also that the bootstrap support values for nodes on the tree must be computed ``manually'', i.e., the program should be repeatedly run on bootstrap samples. This feature is soon planned to be added to the code, but is not currently available.



<br><br>

## Future work and software comparison

The composite likelihood calculation implemented here matches that in PhyNEST (Kong et al., Syst. Biol.  74(1): 53–69, 2025) when the number of reticulations is set to 0 and the data do not include gaps. It does not match that implemented in qAge in PAUP* because the PAUP* implementation further reduced the site pattern counts down to 9 or 11 categories to increase computational efficiency on fixed trees (see paper referenced above for details).

Several extensions of this method are under active development in the Kubatko Research Group, with planned incorporation in this software. In particular, methods for testing hypotheses using composite likelihood ratio tests are being developed and a Bayesian species tree inference method using composite likelihood is in the works. In addition, the following more minor additions to the software are planned:

<ul>
<li> Remove requirement for the tree in the file treefile.tre to have branch lengths</li>
<li> Bootstrap support values on nodes
<li> Additional options for handling gaps and ambiguity codes (currently ignored at the quartet level)
<li> Estimation of the effective population size parameter $\theta$
<li> Posterior site pattern probabilities (Richards and Kubatko, 2021)
<li> Site pattern probabilities under a relaxed clock (Richards and Kubatko, 2022)
<li> Site pattern probabilities for gene tree estimation
</ul>
