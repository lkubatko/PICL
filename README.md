# PICL: Phylogenetic Inference with Composite Likelihood
<b>A package for inferring phylogenies under the MSC and other models</b>

The PICL package includes a collection of methods for analyzing sequence data using composite likelihood. It includes methods for estimating species trees under the multispecies coalescent (MSC). For analyses 
under the MSC, both multilocus and SNP data can be analyzed, and rate variation among sites can be incorporated. Multiple individuals per species can be included. Details about the method can be found in the following publication:

Kubatko, L., S. Kong, E. Webb, and Z. Chen. 2025. The promise of composite likelihood for species-level phylogenomic inference, Evolutionary Journal of the Linnean Society [<a href="https://academic.oup.com/evolinnean/advance-article/doi/10.1093/evolinnean/kzaf008/8127126?login=true">web link</a>]




## Installation

PICL is written in C. After cloning the repository or downloading the code, PICL will need to be compiled.

On a Mac, open a terminal window, navigate to the directory where the source code is stored, and issue the command: gcc main.c -lm -o picl. The program can then be called with `./picl’.

It is recommended that you run the two examples below before trying your own data. The first example will take about 30 minutes to run (reducing the number of optimization steps will reduce the time). The second should run in under 5 minutes.

<br>
<b> Please see the Documentation for additional details.</b>
