# Scavenger biodiversity metrics driving carrion consumption
This repository contains data and R code to reproduce the statistical analyses showed in the paper "Large home range scavengers support higher rates of carcass removal". 

Using five study areas in Spain and South Africa, characterised by different levels of abundance and richness of vultures, large carnivores and other facultative scavenger species, we aim to identify the scavenger community attributes and ecological mechanisms that drive carrion consumption.

More specifically, we test three competitive hypotheses to assess if carrion consumption in vertebrate scavenger communities depends on: i) the presence of key dominant traits (functional identity hypothesis), ii) functional diversity that promotes niche complementarity (functional diversity hypothesis), or iii) the accumulation of individuals and species, irrespective of their trait representation (functional equivalence).

## R files description:

* **0_FD_functions.R**: Functions to estimate Functional Diversity (FD) metrics
* **0_quality_funct_space.R**: R function for computing the quality of functional dendrogramm and multidimensional functional spaces. This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)
* **1_biometric estimation.R**: Code to estimate scavenger biodiversity metrics (functional identity, functional diversity, functional equivalence).
* **2_descriptive_plots.R**: Exploratory correlations of scavenger biodiversity metrics and carrion consumption within each area
* **3_scav_models.R**: GLS models exploring the relationship between scavenger biodiversity metrics and carrion consumption for all areas

## Original data
* **traits.txt**: functional traits for 36 scavenger species
* **dat_final.txt**: carrion consumption rates and scavenger biodiversity metrics for each study site (n=149)


## Original article:

Please, use this citation to reference the data and R code:

```
Gutiérrez-Cánovas, C., Moleón, M., Mateo-Tomás, P, Olea, P.P., Sebastián-González, E., 
Sánchez-Zapata, J.A. Accepted. Large home range scavengers support higher rates of 
carcass removal. Functional Ecology.
```

Please, send questions or problems related with the use of this code to Cayetano Gutiérrez Cánovas (cayeguti@um.es).
