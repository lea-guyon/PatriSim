#!/bin/bash

STARTTIME=$(date +%s)

#############################
######## PARAMETERS #########
#############################

dir=~ # path
burnin=true  # run the burnin (need to be run only once)

chr_size="c(1e6, 1e6, 1e4)"
random_fission=false # true or false
transmission="full" # "full" or "half"
fission_threshold=150
pM=0  # probability for a group to move to another village after a split
violence=false # true or false
descent="unilineal" # "unilineal" or "bilateral"
descent_rule="patrilineal" # "patrilineal" or "matrilineal"
residence_rule_1="multilocal" # "patrilocal" or "matrilocal"
residence_rule_2="patrilocal" # "patrilocal" or "matrilocal"
nb_villages=5
nb_groups=1 # nb of descent groups -> does not make sense for bilateral descent but usefull to normalize villages' sizes
K=300 # carrying capacity per group
declare -i K_total=$nb_villages*$nb_lin*$K # total carrying capacity

########## EXTINCTION RATE IN CASE OF VIOLENCE ###########
e=0.15
##########################################################

sigma=0.1 # variance of the normal law used to draw growth rates
growth_rate=0.01 # growth rate of villages and outgroup, if 0 : population has a constant size
sample_size=20
nbsimu=200 # nb of simulations
cores=40
nameDir="multilocal2patrilineal" # name of the output directory

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
	path=$descent/extended/r=$growth_rate
	rm -rf $dir/Tables/metrics/$path/$nameDir
	mkdir -p $dir/Tables/metrics/$path/$nameDir
	echo "'Replicat'	'Generation'	'mean_nb_children_per_couple'	'var_nb_children_per_couple'	'fathers'	'mothers' 'singleInd'" > $dir/Tables/metrics/$path/$nameDir/metrics.txt
else
	if [ $vl = "T" ]; then
		path=$descent/extended/r=$growth_rate/sigma=$sigma/FT=$fission_threshold/e=$e
	else
		path=$descent/extended/r=$growth_rate/sigma=$sigma/FT=$fission_threshold
	fi
	rm -rf $dir/Tables/metrics/$path/$nameDir
	mkdir -p $dir/Tables/metrics/$path/$nameDir
	echo "'Replicat'    'Generation'    'Nb_of_fissions'    'Nb_of_extinctions'    'Nb_of_lineages'    'fathers'    'Nb_indiv_per_lin'    'var_nb_ind_per_lin'    'Nb_women_per_lin'    'mothers'    'failed_couples'    'Nb_children_per_couple'    'var_nb_children_per_couple'    'mean_lin_depth'    'var_lin_depth'    'mean_migrant_ratio'    'var_migrant_ratio'" > $dir/Tables/metrics/$path/$nameDir/metrics.txt
	echo "'Replicat'	'Generation'	'Village'	'Lineage'	'Nb_children'	'Growth_rate'	'Proba_lin'	'Size'" > $dir/Tables/metrics/$path/$nameDir/nChildrenPerCouple.txt
fi

cd simulations/$path/$nameDir/

## Replace parameters in the slim file ##

cat $dir/SLiM_scripts/bilateral_2_patrilineal.slim | sed "s/bash_Num_villages/${nb_villages}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_random_fission/${rf}/g;s/bash_fission_threshold/${fission_threshold}/g;s/bash_pM/${pM}/g;s/bash_violence/${vl}/g;s/bash_extinction_rate/${e}/g;s/bash_descent_rule/${descent_rule}/g;s/bash_residence_rule_1/${residence_rule_1}/g;s/bash_residence_rule_2/${residence_rule_2}/g;s/bash_sigma/${sigma}/g;s/bash_growth_rate/${growth_rate}/g;s/bash_transmission/${transmission}/g;s/bash_namedir/${nameDir}/g" > "islandmodel.slim"

## Create a new file for each simulation ##
for i in $(seq 1 1 $nbsimu)
do	
	mkdir $dir/simulations/$path/$nameDir/$i
done

cd $dir/simulations/$path/$nameDir/
if $burnin; then
	echo "burnin"
	for i in $(seq 1 1 $nbsimu)
	do	
		cd $dir/simulations/$path/$nameDir/$i
		cat $dir/SLiM_scripts/burnin.slim | sed "s/bash_Num_villages/${nb_villages}/g;s/bash_Num_lin_per_v/${nb_lin}/g;s/bash_total_carrying_capacity/${K_total}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_descent/${descent}/g;s/bash_Num_replicat/${i}/g" > "burnin_${i}.slim"
		echo "slim $i/burnin_${i}.slim"
		cd ..
	done > launcher.txt
	parallel -a launcher.txt -j $cores
fi

echo "SLiM simulations"
cd $dir/simulations/$path/$nameDir

for i in $(seq 1 1 $nbsimu)
do
	cd $dir/simulations/$path/$nameDir/$i
	cat ../islandmodel.slim | sed "s/bash_Num_replicat/${i}/g" > "islandmodel_${i}.slim"
    echo "slim $i/islandmodel_${i}.slim > $i/outputSlim${i}.slim"
		
	cd ..
done > launcher.txt
parallel -a launcher.txt -j $cores > outputErrors.slim

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
	echo "python $dir/Python_scripts/subset_trees_villages_bilateral_descent.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size $sample_size -K $K_total -o $dir/simulations/$path/$nameDir/$i/ > $i/outputPy${i}.txt"
	echo "python $dir/Python_scripts/subset_trees_villages_unilineal_descent.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size $sample_size -K $K_total -d $descent_rule -o $dir/simulations/$path/$nameDir/$i/ -t $dir/Tables/metrics/$path/$nameDir > $i/outputPy${i}.txt"
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
echo "It takes $(($ENDTIME - $STARTTIME3)) seconds to complete this task"

##################################
####### BAYESIAN INFERENCE #######
##################################

echo "write .nex files"
STARTTIME3=$(date +%s)

cd $dir/BEAST/$descent/extended/
rm -rf $dir/BEAST/$path/$nameDir
mkdir -p $dir/BEAST/$path/$nameDir

for i in $(seq 1 1 $nbsimu)
do	
	mkdir $dir/BEAST/$path/$nameDir/$i
done

echo "Generate .nex files"
for i in $(seq 1 1 $nbsimu)
do
	echo "python $dir/Python_scripts/write_nexus.py -s $dir/simulations/$path/$nameDir/ -rep $i --sample-size 20 -K $K_total -gen 200 -o $dir/BEAST/$path/$nameDir/$i/ > $dir/simulations/$path/$nameDir/$i/outputPyNex${i}.txt"
done > launcher.txt
parallel -a launcher.txt -j $cores

cd $dir/BEAST/
# transform nexus
echo "Transform .nex files"
for i in $(seq 1 1 $nbsimu)
do
	echo "bash transform_nexus.txt -d $dir -p $path -n $nameDir -R ${i} -g 200"
done > launcher.txt
parallel -a launcher.txt -j $cores

cd $dir/BEAST/
# create xml file
echo "Generate .xml files"
for i in $(seq 1 1 $nbsimu)
do
	echo "beastgen template_Mito.xml $path/$nameDir/$i/Sim_Mito.nex $path/$nameDir/$i/Sim_Mito.xml"
	echo "beastgen template_Y.xml $path/$nameDir/$i/Sim_Y.nex $path/$nameDir/$i/Sim_Y.xml"
done > launcher.txt
parallel -a launcher.txt -j $cores

cd $dir/BEAST/$path/$nameDir/
# beast inference
echo "bayesian inference"
for i in $(seq 1 1 $nbsimu)
do
	echo "beast -working -beagle_scaling always -overwrite $i/Sim_Mito.xml"
	echo "beast -working -beagle_scaling always -overwrite $i/Sim_Y.xml"
done > launcher.txt
parallel -a launcher.txt -j $cores

echo "Generate Skyline Plots"
for i in $(seq 1 1 $nbsimu)
do
	echo "Rscript $dir/R_scripts/skyline_plots.R '$dir/BEAST/$path/$nameDir/$i/'"
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME3)) seconds to complete this task"