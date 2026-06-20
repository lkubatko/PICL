# PICL: Phylogenetic Inference with Composite Likelihood
<b>A package for inferring phylogenies under the MSC and other models</b>

The PICL package includes a collection of methods for analyzing sequence data using composite likelihood. It includes methods for estimating species trees under the multispecies coalescent (MSC). For analyses 
under the MSC, both multilocus and SNP data can be analyzed, and rate variation among sites can be incorporated. Multiple individuals per species can be included. Details about the method can be found in the following publication:

Kubatko, L., S. Kong, E. Webb, and Z. Chen. 2025. The promise of composite likelihood for species-level phylogenomic inference, Evolutionary Journal of the Linnean Society [<a href="https://academic.oup.com/evolinnean/advance-article/doi/10.1093/evolinnean/kzaf008/8127126?login=true">web link</a>]

<br>

## jPICL -- a GUI for the PICL program
<a href="https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/people/daniel-huson/">Daniel Huson</a> has created a GUI for the PICL program, called jPICL. You can download jPICL <a href="https://github.com/husonlab/jpicl/releases">here.</a>  We recommend getting started with PICL using jPICL. Shorter runs
can be carried out entirely within the GUI. For longer runs or use in an HPC environment, jPICL can be used to generate a `settings` file that is compatible with PICL.

<br>

## Installation from source code

If you would like to compile PICL directly on your system, you can clone the repository or download the code. PICL will then need to be compiled.   The MacOS version of PICL uses routines from the <a href="https://github.com/boostorg/boost">Boost C++ libraries</a>. These will need to installed prior to compilation (see the Documentatino for instructions).

After downloading the soruce code, open a terminal window, navigate to the directory where the source code is stored, and issue the command: `make`. The program can then be called with `./picl`.

<br><br>

<b> Please see the Documentation for additional details.</b>
