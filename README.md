# The evolution of inversions suppressing recombination between sex chromosomes

## Overview

This is a GitHub repository for the development of a theoretical population genetics research project that is now published under the title "*The evolution of suppressed recombination between sex chromosomes by chromosomal inversions*" (doi: [https://doi.org/10.1101/2020.03.23.003558](https://doi.org/10.1101/2020.03.23.003558)). Here you can find all of the necessary code to reproduce the simulations presented in the published paper and appendices, and the LaTeX files used to compile the manuscript. Other supplementary material for the paper can be downloaded as an online appendix from [link to publisher website](URL).


## Abstract

The idea that sex-differences in selection drive the evolution of suppressed recombination between sex chromosomes is well-developed in population genetics. Yet, despite a now classic body of theory, empirical evidence that sexually antagonistic selection drives the evolution of recombination arrest remains meagre and alternative hypotheses underdeveloped. Here, we investigate whether the length of 'evolutionary strata' formed by chromosomal inversions expanding the non-recombining sex-linked region (SLR) on proto sex chromosomes can be informative of how selection influenced their fixation. We develop population genetic models to show how the length of an SLR-expanding inversion, and the presence of partially recessive deleterious mutational variation, affect the fixation probability of three different classes of inversions: (*i*) neutral, (*ii*) directly beneficial (i.e., due to breakpoint or positional effects), and (*iii*) those capturing sexually antagonistic (SA) loci. Our models indicate that neutral inversions, and those capturing an SA locus in linkage disequilibrium with the ancestral SLR, will exhibit a strong fixation bias towards small inversions; while unconditionally beneficial inversions, and those capturing a previously unlinked SA locus, will generally favour fixation of larger inversions. The signatures of evolutionary strata size left behind by different selection regimes are strongly influenced by parameters affecting the deleterious mutation load, the physical position of the ancestral SLR, and the shape of the length distribution of new inversions (about which we know very little).

## Citing information

*Preprint citation*:

Olito, C. and J.K. Abbott. 2020. The evolution of suppressed recombination between sex chromosomes by chromosomal inversions. bioRXiv doi: [https://doi.org/10.1101/2020.03.23.003558](https://doi.org/10.1101/2020.03.23.003558).

*Paper citation*:

Citing information for the final paper will be provided when it is made [available through the publisher](URL). You can also contact me directly if you would like a reprint. 

##  Instructions

This repository provides all code necessary to (1) rerun the simulations and (2) produce the main figures from the paper as .pdf's. To do so, please follow these basic steps:

1. Clone the repo using the following: `git clone https://https://github.com/colin-olito/inversionSize-ProtoSexChrom`. Alternatively, on the project main page on GitHub, click on the green button `clone` or `download` and then click on `Download ZIP`.  
2. Check that you have a recent version of [`R`](https://www.r-project.org/) installed. 
3. Make sure that the working directory for your R session is the root directory of this repo (e.g., `inversionSize-ProtoSexChrom-master/`).
4. Run `./R/run-WFSims-prFix-recDel.R` either interactively in R or in terminal. The simulations can take quite some time to run.
5. *Note*: We use CM fonts in the figures. To do this, be sure to correctly install the `R` font packages `extrafont` and `fontcm`. Alternatively, comment out L.4-6 in `./R/functions-figures`, and change the default font family to 'Arial' by swapping L.27 & L.28.
6. Run `makeFigs.R`, which will read the simulation output files and generate the main figures in the paper and supplementary material.  


## Repository structure and contents 

The directories/files in this repository needed to reproduce the results for this study are as follows:  

- **`R`**   
	- `functions-beneficial-recDel.R`  
	- `functions-figures.R`  
	- `functions-SA-recDel.R`  
	- `functions-sheltering-recDel.R`  
- **`output`**   
	- **`data`**   
			- **`simResults`**   
				- `SA-eqFreqs-*.R`  
				- `SA-sI-*.R`  
	- **`figures`***  
- `makeFigs.R`  
- `run-WFSims-prFix-recDel.R`  
- `LICENSE.txt`   

**Note:** * `./output` and `./output/figures` directories will be created locally the first time `run-WFSims-prFix-recDel.R` is run (*if needed*).


### File & variable descriptions

Plotting function files
- `functions-figures.R`: general plotting functions.   

Simulation function files for each selection scenario.
- `functions-beneficial-recDel.R`  
- `functions-figures.R`  
- `functions-SA-recDel.R`  
- `functions-sheltering-recDel.R`  

Executables
- `run-WFSims-prFix-recDel.R`: executable functions for W-F simulations of inversions under deleterious mutation pressure.   
- `makeFigs.R`: executable functions to create .pdf figures using simulation output files.

License    
- `LICENSE.txt`: MIT license for this repository.  


## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the inversionSize github [issues page](https://github.com/colin-olito/inversionSize-ProtoSexChrom/issues), or inform me directly by sending a brief email detailing the problem you encountered to colin.olito at biol dot lu dot se.