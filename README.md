# README

Pipeline simulating human genomes under scenarios with different descent (bilateral, patrilineal, matrilineal) and residence (patrilocal, matrilocal) rules and computing diversity estimators on simulated genetic data.

## Installation
`git clone https://github.com/lea-guyon/PatriSim`

## Summary

To run the pipeline, modify the `PARAMETERS` section in the file `main.sh` (see Fig 1 in article). See "Parameters information" to have a description of each parameter. Then execute the .sh file with bash in a terminal. 

`bash main.sh`

This pipeline:
1) simulates autosomes, X chromosomes, Y chromosomes and mtDNA under the scenario specified by the parameters using the sofware SLiM (see "SLiM models"). 
2) converts SLiM outputs into VCF files with custom python scripts (see "Output VCF")
3) computes nucleotide diversity using custom python scripts (see "Nucleotide diversity")
4) generates BEAST output to draw Bayesian Skyline Plots (see "Bayesian inference")
5) Plot BSPs using Art Poon's R script: "plot-skyline.R" (https://gist.github.com/ArtPoon/4c6cb4bdae02bf4bbd659b131a2afa9e)

NB: if you already ran the script once, it is recommended to set the `burnin` parameter to "false".

To run the scenarios presented in the paper, use the parameters values indicated in `scenarios_parameters.md`

## Requirements

### System requirements
The pipeline requires a standard computer with enough RAM to support the in-memory operations.

#### OS requirements
The pipeline is designed to be run on Linux. It has been tested on Ubuntu 20.04.

### SLiM
The pipeline requires SLiM v4.0.1. To install SLiM, follow the instructions given in the manual (https://messerlab.org/slim/). 

### Python dependencies
To install python (3.8) requirements enter:  
`pip3 install -r requirements.txt`

## Parameters information

`dir`: input path

`chr_size`: vector of chromosome sizes. First element: autosomes and X. Second element: Y chromosome. Third element: mtDNA.

`random_fission`: true or false. If true, the fission type is random.

`transmission`: "full" or "half". If "full", the growth rate of a splitting descent group is transmitted to the two resulting descent groups after the fission event. If "half" it is transmitted to one out of the two new groups.

`fission_threshold`: number of individuals from which a descent group splits

`pM`: float between 0 and 1. pM is the probability for the smallest group resulting from a split to move to another village.

`violence`: true or false. If true, there is violent intergroup competition between descent groups.

`descent`: "unilineal" or "bilateral"

`descent_rule`: "patrilineal" or "matrilineal"

`residence_rule_1`: "patrilocal" or "matrilocal" or "na", corresponding to the first residence rule implemented in the "bilateral2patrilineal" model

`residence_rule_2`: "patrilocal" or "matrilocal" or "na", corresponding to the second residence rule implemented in the "bilateral2patrilineal" model

`nb_villages`: number of villages in the population

`nb_groups`: number of descent groups -> does not make sense for bilateral descent but usefull to normalize villages' size

`K`: initial number of individuals per descent group

`polygyny`: "F" or "T". If "T", males mate with several women (change the slim code to modify the rate of polygyny)

`e`: extinction rate. Proportion of males killed at each generation.

`mf`: proportion of females migrating to a new village at each generation.

`mm`: proportion of males migrating to a new village at each generation.

`sigma`: variance of the normal law used to draw relative fitnesses

`growth_rate`: growth rate of villages and outgroup, if 0 : population has a constant size

`sample_size`: number of individuals to sample in each village

`nbsimu`: number of replicates

`cores`: number of cores

`nameDir`: name of the output directory

## SLiM models
*Contains .slim files modeling populations with different descent (patrilineal, matrilineal, bilateral) and post-marital residence rules (patrilocal, matrilocal)*

- **burnin.slim**: simulate a panmictic population made of `K_total` individuals for 20,000 generations. This preliminary simulation is taken in entry of the models of villages with bilateral or unilineal descent. 

- **bilateral_descent.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. See Fig 1.

- **unilineal_descent.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals divided into `nb_groups` descent groups at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. The descent rule is defined by `descent_rule`. See Fig 1.

- **bilateral2patrilineal.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. Then there is a transition to a patrilineal system. This continues for another 100 generations. See Supplementary Fig. 2.

- **patrilineal2bilateral.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals divided into `nb_groups` descent groups at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. The descent rule is defined by `descent_rule`. Then there is a transition to a bilateral system. This continues for another 100 generations.

- **patrilineal_strict2relaxed_.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals divided into `nb_groups` descent groups at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. The descent rule is defined by `descent_rule`. The parameter `sigma` controlling the variance in reproductive success is set to 0.1. Then there is a transition to a patrilineal system without variance in reproductive success. This continues for another 100 generations.

## Output VCF
*.py scripts taking .trees files (generated by simulating populations in SLiM) in entry to produce .vcf files for each chromosome and each village in every simulation*

In `Python_scripts` folder:
- **subset_trees_villages_bilateral_descent.py**
- **subset_trees_villages_unilineal_descent.py**

## Nucleotide diversity
*.py files used to compute diversity metrics (Pi) from .vcf files*

In `Python_scripts` folder:
- **metrics.py** : file containing functions to compute alleles frequencies and diversity metrics

- **Pi_villages.py** : compute nucleotidic diversity for models of villages

- **Global_Pi_villages.py** : compute global nucleotidic diversity for models of villages

## Bayesian inference
*generates .nex BEAST inputs and run beauti and BEAST softwares*

In `Python_scripts` folder:
- **write_nexus_bilateral.py**
- **write_nexus.py**

In `BEAST` folder:

- **transform_nexus.sh**
- **template_Y.xml**
- **template_Mito.xml**

In `R_scripts` folder:

-**skyline_plot.R**: adapted from Art Poon's github. Draw Bayesian Skyline Plots.
