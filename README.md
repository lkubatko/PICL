# PICL: Phylogenetic Inference with Composite Likelihood
<b>A package for inferring phylogenies under the MSC and other models</b>

The PICL package includes a collection of methods for analyzing sequence data using composite likelihood. It includes methods for estimating species trees under the multispecies coalescent (MSC). For analyses 
under the MSC, both multilocus and SNP data can be analyzed, and rate variation among sites can be incorporated. Multiple individuals per species can be included. Details about the method can be found in the following publication:

Kubatko, L., S. Kong, E. Webb, and Z. Chen. 2025. The promise of composite likelihood for species-level phylogenomic inference, Evolutionary Journal of the Linnean Society [<a href="https://academic.oup.com/evolinnean/advance-article/doi/10.1093/evolinnean/kzaf008/8127126?login=true">web link</a>]


## jPICL -- a GUI for the PICL program
<a href="https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/people/daniel-huson/">Daniel Huson</a> has created a GUI for the PICL program, called jPICL. You can download jPICL <a href="https://github.com/husonlab/jpicl/releases">here.</a>  

## Installation from source code

PICL is written in C. After cloning the repository or downloading the code, PICL will need to be compiled.   The MacOS version of PICL uses routines from the <a href="https://github.com/boostorg/boost">Boost C++ libraries</a>. These will need to installed prior to compilation (see below for instructions).

Open a terminal window, navigate to the directory where the source code is stored, and issue the command: `make`. The program can then be called with `./picl`.

<br><br>

<b> Please see the Documentation for additional details.</b>
