To reproduce the results obtained in the paper, run the script `main.sh` with the following parameters:

## Scenario 1
chr_size="c(1e6, 1e6, 1e4)"  
descent="bilateral"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  \
mf=0.1 \
mm=0  
growth_rate=0.01   
sample_size=20  
nbsimu=200  
cores=30  
nameDir="patrilocal_villages"

## Scenario 2a
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=true  
transmission="full"  
fission_threshold=150  \
pM=0  \
friendlyFission="T" \
violence=false  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F"  \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0 \
sigma=0  
growth_rate=0.01  
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_a"  

## Scenario 2b
chr_size="c(1e6, 1e6, 1e4)"   
random_fission=true   
transmission="full"  
fission_threshold=150 \
pM=0 \
friendlyFission="T"  
violence=true  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0  
sigma=0  
growth_rate=0.01   
sample_size=20   
nbsimu=200  
cores=40  
nameDir="patrilineal_villages_b"  

## Scenario 2c
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150 \
pM=0 \
friendlyFission="T"  
violence=false  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K  

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1
mm=0  
sigma=0  
growth_rate=0.01  
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_c"  

## Scenario 2d
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150 \
pM=0 \
friendlyFission="T"  
violence=true  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0 \
sigma=0 \
growth_rate=0.01  
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_e"

## Scenario 2e
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=true  
transmission="full"  
fission_threshold=150 \
pM=0 \
friendlyFission="T"  
violence=false  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0 \
sigma=0.1  
growth_rate=0.01  
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_d"  

## Scenario 2f
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=true  
transmission="full"  
fission_threshold=150 \
pM=1 \
friendlyFission="T"  
violence=true  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0 \
sigma=0.1  
growth_rate=0.01   
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_f"  

## Scenario 2g
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150 \
pM=0 \
friendlyFission="T"  
violence=false  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0 \
sigma=0.1  
growth_rate=0.01  
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_g"  

## Scenario 2h
chr_size="c(1e6, 1e6, 1e4)"  
random_fission=false  
transmission="full"  
fission_threshold=150 \
pM=0 \
friendlyFission="T"  
violence=true  
descent="unilineal"  
descent_rule="patrilineal"  
nb_villages=5  
nb_groups=3  
K=100  
polygyny="F" \
declare -i K_total=$nb_villages*$nb_groups*$K   

######## EXTINCTION RATE IN CASE OF VIOLENCE ########  \
e=0.15   
#############################################

mf=0.1 \
mm=0 \
sigma=0.1  
growth_rate=0.01  
sample_size=20  
nbsimu=200   
cores=40  
nameDir="patrilineal_villages_h"