#!/bin/bash

STARTTIME=$(date +%s)

#############################
######## PARAMETERS #########
#############################

dir=~ # path

chr_size="c(1e6, 1e6, 1e4)"
random_fission=false # true or false
transmission="full" # "full" or "half"
fission_threshold=100
friendlyFission="T" # "T" or "F"
violence=false # true or false
descent="unilineal" # "unilineal" or "bilateral"
descent_rule="patrilineal" # "patrilineal" or "matrilineal"
residence_rule="patrilocal" # "patrilocal" or "matrilocal"
nb_villages=5
nb_groups=3 # nb of descent groups -> does not make sense for bilateral descent but usefull to normalize villages' sizes
K=100 # carrying capacity per lineage
polygyny="F" # "F" or "T"
declare -i K_total=$nb_villages*$nb_groups*$K # total carrying capacity

########## EXTINCTION RATE IN CASE OF VIOLENCE ###########
e=0.15
##########################################################

migration_rate=0.05
sigma=0.1 # variance of the normal law used to draw growth rates
growth_rate=0.01 # growth rate of villages and outgroup, if 0 : population has a constant size
sample_size=20
nbsimu=200 # nb of simulations
cores=40
nameDir="patrilineal_villages" # name of the output directory

############################
####### SIMULATIONS ########
############################

if $random_fission; then
	rf="T"
else
    rf="F"
fi
if $violence; then
	vl='T'
else
    vl="F"
fi

echo "Starting simulations"

cd $dir

if [ "$descent" = "bilateral" ]; then
	path=$descent/regular/r=$growth_rate
	rm -rf $dir/Tables/metrics/$path/$nameDir
	mkdir -p $dir/Tables/metrics/$path/$nameDir
	echo "'Replicat'	'Generation'	'mean_nb_children_per_couple'	'var_nb_children_per_couple'	'fathers'	'mothers' 'singleInd'" > $dir/Tables/metrics/$path/$nameDir/metrics.txt
else
	if [ $vl = "T" ]; then
		path=$descent/regular/r=$growth_rate/sigma=$sigma/FT=$fission_threshold/e=$e
	else
		path=$descent/regular/r=$growth_rate/sigma=$sigma/FT=$fission_threshold
	fi
	rm -rf $dir/Tables/metrics/$path/$nameDir
	mkdir -p $dir/Tables/metrics/$path/$nameDir
	echo "'Replicat'    'Generation'    'N_ind'    'Nb_of_fissions'    'Nb_of_extinctions'    'Nb_of_lineages'    'fathers'    'Nb_indiv_per_lin'    'var_nb_ind_per_lin'    'Nb_women_per_lin'    'mothers'    'failed_couples'    'singleInd'	'Nb_children_per_couple'    'var_nb_children_per_couple'    'mean_lin_depth'    'var_lin_depth'    'mean_migrant_ratio'    'var_migrant_ratio' 'meanFissionTime' 'varFissionTime'" > $dir/Tables/metrics/$path/$nameDir/metrics.txt
	echo "'Replicat'	'Generation'	'Lineage'	'nChildren'" > $dir/Tables/metrics/$path/$nameDir/nChildrenPerCouple.txt
	echo "'Replicat'    'Generation'    'GroupDepth'" > $dir/Tables/metrics/$path/$nameDir/groupDepth.txt
	echo "'Replicat'    'Generation'    'FissionTime'" > $dir/Tables/metrics/$path/$nameDir/fissionTime.txt
fi

cd simulations/$path/$nameDir/

## Replace parameters in the slim file ##

if [ "$descent" = "bilateral" ]; then
	cat $dir/SLiM_models/bilateral_descent.slim | sed "s/bash_Num_villages/${nb_villages}/g;s/bash_chr_size/${chr_size}/g;s/bash_migration_rate/${migration_rate}/g;s/bash_residence/${residence_rule}/g;s/bash_growth_rate/${growth_rate}/g;s/bash_namedir/${nameDir}/g;s/bash_polygyny/${polygyny}/g" > "islandmodel.slim"
else
    cat $dir/SLiM_models/unilineal_descent.slim | sed "s/bash_Num_villages/${nb_villages}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_random_fission/${rf}/g;s/bash_fission_threshold/${fission_threshold}/g;s/bash_friendly_fission/${friendlyFission}/g;s/bash_violence/${vl}/g;s/bash_extinction_rate/${e}/g;s/bash_migration_rate/${migration_rate}/g;s/bash_descent_rule/${descent_rule}/g;s/bash_residence_rule/${residence_rule}/g;s/bash_sigma/${sigma}/g;s/bash_growth_rate/${growth_rate}/g;s/bash_transmission/${transmission}/g;s/bash_namedir/${nameDir}/g" > "islandmodel.slim"
fi

## Create a new file for each simulation ##
echo "burnin"
for i in $(seq 1 1 $nbsimu)
do	
	mkdir $dir/simulations/$path/$nameDir/$i
	cd $dir/simulations/$path/$nameDir/$i

	# comment to skip burn-in simulation
	cat $dir/SLiM_models/burnin.slim | sed "s/bash_Num_villages/${nb_villages}/g;s/bash_nGroupsPerVillage/${nb_groups}/g;s/bash_total_carrying_capacity/${K_total}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_descent/${descent}/g;s/bash_Num_replicat/${i}/g" > "burnin_${i}.slim"
	echo "slim $i/burnin_${i}.slim"
	cd ..
done > launcher.txt
parallel -a launcher.txt -j $cores

## Run simulations ##
echo "run simulations"
cd $dir/simulations/$path/$nameDir

for i in $(seq 1 1 $nbsimu)
do
	cd $dir/simulations/$path/$nameDir/$i
	cat ../islandmodel.slim | sed "s/bash_Num_replicat/${i}/g" > "islandmodel_${i}.slim"

    echo "slim $i/islandmodel_${i}.slim > $i/outputSlim${i}.slim"
		
	cd ..
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

##################################
####### OUTPOUT VCF FILES ########
##################################

cd $dir/simulations/$path/$nameDir/
echo "output VCF files"
STARTTIME2=$(date +%s)

for i in $(seq 1 1 $nbsimu)
do
	cd $dir/simulations/$path/$nameDir/

	if [ "$descent" = "bilateral" ]; then
        echo "python $dir/Python_scripts/subset_trees_villages_bilateral_descent.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size $sample_size -K $K_total -o $dir/simulations/$path/$nameDir/$i/ > $i/outputPy${i}.txt"
	else
		echo "python $dir/Python_scripts/subset_trees_villages_unilineal_descent.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size $sample_size -K $K_total -d $descent_rule -o $dir/simulations/$path/$nameDir/$i/ -t $dir/Tables/metrics/$path/$nameDir > $i/outputPy${i}.txt"
	fi
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME2)) seconds to complete this task"

##################################
######## COMPUTE METRICS #########
##################################

echo "compute diversity metrics"
STARTTIME4=$(date +%s)
cd $dir/

rm -rf $dir/Tables/Pi/$path/$nameDir
mkdir -p $dir/Tables/Pi/$path/$nameDir
python Python_scripts/Pi_villages.py -s $dir/simulations/$path/$nameDir/ -o $dir/Tables/Pi/$path/$nameDir/

python Python_scripts/Global_Pi_villages.py -s $dir/simulations/$path/$nameDir/ -o $dir/Tables/Pi/$path/$nameDir/

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME4)) seconds to complete this task"

##################################
####### BAYESIAN INFERENCE #######
##################################

echo "write .nex files"
STARTTIME3=$(date +%s)

cd $dir/BEAST/$descent/regular/
rm -rf $dir/BEAST/$descent/regular/r=$growth_rate/$nameDir
mkdir -p $dir/BEAST/$descent/regular/r=$growth_rate/$nameDir

for i in $(seq 1 1 $nbsimu)
do	
	mkdir $dir/BEAST/$descent/regular/r=$growth_rate/$nameDir/$i
done


for i in $(seq 1 1 $nbsimu)
do
	if [ "$descent" = "bilateral" ]; then
        echo "python $dir/Python_scripts/write_nexus_bilateral.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size 50 -K $K_total -gen 100 -o $dir/BEAST/$descent/regular/r=$growth_rate/$nameDir/$i/ -t $dir/Tables/metrics/$path/$nameDir/ > $dir/simulations/$path/$nameDir/$i/outputPyNex${i}.txt"
	else
		echo "python $dir/Python_scripts/write_nexus.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size 50 -K $K_total -gen 100 -o $dir/BEAST/$descent/regular/r=$growth_rate/$nameDir/$i/ > $dir/simulations/$path/$nameDir/$i/outputPyNex${i}.txt"
	fi
done > launcher.txt
parallel -a launcher.txt -j $cores

cd $dir/BEAST/
# transform nexus
for i in $(seq 1 1 $nbsimu)
do
	echo "bash transform_nexus.txt -d $dir -D $descent -t regular -r ${growth_rate} -n $nameDir -R ${i} -g 100"
done > launcher.txt
parallel -a launcher.txt -j $cores

# create xml file
for i in $(seq 1 1 $nbsimu)
do
	echo "beastgen template_Mito.xml $descent/regular/r=$growth_rate/$nameDir/$i/Sim_Mito.nex $descent/regular/r=$growth_rate/$nameDir/$i/Sim_Mito.xml"
	echo "beastgen template_Y.xml $descent/regular/r=$growth_rate/$nameDir/$i/Sim_Y.nex $descent/regular/r=$growth_rate/$nameDir/$i/Sim_Y.xml"
done > launcher.txt
parallel -a launcher.txt -j $cores

cd $dir/BEAST/$descent/regular/r=$growth_rate/$nameDir/

# beast inference
echo "bayesian inference"
for i in $(seq 1 1 $nbsimu)
do
	echo "beast -working -beagle_scaling always -overwrite $i/Sim_Mito.xml"
	echo "beast -working -beagle_scaling always -overwrite $i/Sim_Y.xml"
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME3)) seconds to complete this task"