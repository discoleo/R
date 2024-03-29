

# Modeling Infectious Diseases


## Introduction

**TODO: ...**

## Work

### Epidemic Simulator
* based on the code developed for the BSc thesis (Anamaria Bica);
* improved code and redesigned as R package;
* GitHub: https://github.com/discoleo/EpidemicSimulator

### Bachelor Thesis: Anamaria Bica (co-supervisor: L. Mada);
* partly based on an initial group-project for students (initiated/supervised by L. Mada);
* the initial code was redesigned, standardized and extended with various new compartment models;
* GitHub: [old] https://github.com/BicaAnamaria/EpidemicSimulator



## R Packages

### Epidemiology

#### Covid / Coronavirus:
* COVID19, covidcast, covid19.analytics, covidregionaldata, coronavirus, nCov2019, covid19dbcand, canadacovid, covid19us, covid19india, covid19srilanka, covid19br, COVIDIBGE (Brazil), covid19sf (San Francisco), oxcgrt (Gov responses);
* Europe: covid19france, covid19swiss, UKB.COVID19, covid19italy, covidsymptom (Sweden);

#### Literature:
* pubmed.mineR;
#### Simulations/Modeling:
* EpiModel, PandemicLP, babsim.hospital, cif (for ICU), niaidMI, weibull4, covidprobability, bets.covid19;
* modelSSE: Modelling Infectious Disease Superspreading from Contact Tracing
#### Genetics/Mutations:
* CovidMutations, UKB.COVID19, VERSO (phylo, Bioconductor);
#### Other:
* opitools (opinions on Twitter/YT), SIPDIBGE (Brazil microdata);

Note: Bioconductor packages not yet searched;


### Solvers

1. Package deSolve
Vignette: Solving Initial Value Differential Equations in R.
> https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf

2. Package FME:
Vignette: Inverse Modelling, Sensitivity and Monte Carlo Analysis in R Using Package FME
> https://cran.r-project.org/web/packages/FME/vignettes/FME.pdf

3. Package mosaicCalc
Vignette: Instructors’ Guide
> https://cran.r-project.org/web/packages/mosaicCalc/vignettes/Instructors.html

4. Package deBInfer
Vignette: Bayesian inference for a population growth model of the chytrid fungus.
> https://cran.r-project.org/web/packages/deBInfer/vignettes/chytrid_dede_example.pdf

5. Package rodeo
Vignette: Basic Use and Sample Applications.
> https://cran.r-project.org/web/packages/rodeo/vignettes/rodeoVignette.pdf

6. Package deBif
Vignette: A package for bifurcation analysis of ordinary differential equation systems.
> https://cran.r-project.org/web/packages/deBif/vignettes/deBif-manual.pdf


## Literature

1. Cummings, Derek A.T., et al. Traveling waves in the occurrence of dengue haemorrhagic fever in Thailand.
Nature, vol. 427, no. 6972, 22 Jan. 2004, pp. 344-347.

2. Huang, N. E. et al. The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis.
Proc. R. Soc. Lond. A 454, 903–995 (1998).

3. Balabdaoui, Fadoua & Mohr, Dirk. Age-stratified discrete compartment model of the COVID 19 epidemic with application to Switzerland.
Nature (2020).

4. Romano, S, Fierro, A, Liccardo, A. Beyond the peak: A deterministic compartment model for exploring the Covid-19 evolution in Italy.
Plos One (2020).


### Contact Tracing

1. Johannes Müller. Contact Tracing & Super-Spreaders in the Branching-Process Model.
https://www.youtube.com/watch?v=MCkAahJ08GM

