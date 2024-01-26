#!/bin/bash

STARTTIME=$(date +%s)

#############################
######## PARAMETERS #########
#############################

dir="path"

burnin=false  # run the burnin (need to be run only once)

chr_size="c(1e6, 1e6, 1e4)"
random_fission=false
transmission="full"
fission_threshold=150
pM=0  # probability for a group to move to another village after a split
violence=false
descent="unilineal" # "unilineal" or "bilateral"
descent_rule="patrilineal" # "patrilineal" or "matrilineal"
nb_villages=5
nb_group=3 # nb of groups -> doesn't make sense for bilateral descent but usefull to normalize villages' sizes
K=100 # carrying capacity per group
polygyny="F"
supermale="F"
tPOLY="F"
declare -i K_total=$nb_villages*$nb_group*$K # total carrying capacity
declare -i K_total_o=2*K_total # total carrying capacity in simulations with outgroup

########## EXTINCTION RATE IN CASE OF VIOLENCE ###########
e=0.15
##########################################################

mf=0.1  # female migration rate
mm=0  # male migration rate
sigma=0.1
growth_rate=0.01 #0.01 # growth rate of villages and outgroup, if 0 : population has a constant size
sample_size=20
pseudohap=False
nbsimu=200 # nb of simulations
cores=40

############################
####### SIMULATIONS ########
############################

echo "Starting simulations"

cd $dir
vl="F"
rf="F"
nameDir="patrilineal_strict2relaxed"
path=$descent/extended/r=$growth_rate/sigma=$sigma/FT=$fission_threshold

rm -rf $dir/Tables/Simulations_metrics/$path/$nameDir
mkdir -p $dir/Tables/Simulations_metrics/$path/$nameDir
echo "'Replicat'    'Generation'    'N_ind'    'Nb_of_fissions'    'Nb_of_extinctions'    'Nb_of_groups'    'fathers'    'Nb_indiv_per_group'    'var_nb_ind_per_group'    'Nb_women_per_group'    'mothers'    'failed_couples'    'singleInd'	'Nb_children_per_couple'    'var_nb_children_per_couple'    'mean_group_depth'    'var_group_depth'    'mean_migrant_ratio'    'var_migrant_ratio' 'meanFissionTime' 'varFissionTime'" > $dir/Tables/Simulations_metrics/$path/$nameDir/metrics.txt
echo "'Replicat'	'Generation'	'Village'	'group'	'Nb_children'	'Growth_rate'	'Proba_group'	'Size'" > $dir/Tables/Simulations_metrics/$path/$nameDir/proba_group.txt
echo "'Replicat'	'Generation'	'Village'	'growth_rate'	'var_growth_rate'" > $dir/Tables/Simulations_metrics/$path/$nameDir/growthRates.txt
echo "'Replicat'	'Generation'	'group'	'nChildren'" > $dir/Tables/Simulations_metrics/$path/$nameDir/nChildrenPerCouple.txt
echo "'Replicat'    'Generation'    'GroupDepth'" > $dir/Tables/Simulations_metrics/$path/$nameDir/groupDepth.txt
echo "'Replicat'	'Generation'	'Step'	'nMales'	'nInds'	'sexRatio'" > $dir/Tables/Simulations_metrics/$path/$nameDir/sexRatio.txt

rm -rf Simulations_folders/$path/$nameDir
mkdir -p Simulations_folders/$path/$nameDir

cd $dir/Simulations_folders/$path/$nameDir/

## Replace parameters in the slim file ##

cat $dir/SLiM_scripts/strict2relaxed.slim | sed "s/bash_wd/${dir}/g;s/bash_Num_villages/${nb_villages}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_random_fission/${rf}/g;s/bash_fission_threshold/${fission_threshold}/g;s/bash_pM/${pM}/g;s/bash_violence/${vl}/g;s/bash_extinction_rate/${e}/g;s/bash_mf_ratio/${mf}/g;s/bash_mm_ratio/${mm}/g;s/bash_descent_rule/${descent_rule}/g;s/bash_sd/${sd}/g;s/bash_growth_rate/${growth_rate}/g;s/bash_polygyny/${polygyny}/g;s/bash_t_polygyny/${tPOLY}/g;s/bash_supermale/${supermale}/g;s/bash_transmission/${transmission}/g;s/bash_dir_table/${nameDir}/g" > "islandmodel.slim"


for i in $(seq 1 1 $nbsimu)
do	
	mkdir $dir/Simulations_folders/$path/$nameDir/$i
done

cd $dir/Simulations_folders/$path/$nameDir/
if $burnin; then
	echo "burnin"
	for i in $(seq 1 1 $nbsimu)
	do	
		cd $dir/Simulations_folders/$path/$nameDir/$i
		cat $dir/SLiM_scripts/burnin.slim | sed "s/bash_wd/${dir}/g;s/bash_Num_villages/${nb_villages}/g;s/bash_Num_group_per_v/${nb_group}/g;s/bash_total_carrying_capacity/${K_total}/g;s/bash_carrying_capacity/${K}/g;s/bash_chr_size/${chr_size}/g;s/bash_descent/${descent}/g;s/bash_Num_replicat/${i}/g" > "burnin_${i}.slim"
		echo "slim $i/burnin_${i}.slim"
		cd ..
	done > launcher.txt
	parallel -a launcher.txt -j $cores
fi

echo "SLiM simulations"
cd $dir/Simulations_folders/$path/$nameDir
for i in $(seq 1 1 $nbsimu)
do
	cd $dir/Simulations_folders/$path/$nameDir/$i
	cat ../islandmodel.slim | sed "s/bash_Num_replicat/${i}/g" > "islandmodel_${i}.slim"

    echo "slim $i/islandmodel_${i}.slim > $i/outputSlim${i}.slim"
		
	cd ..
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

##################################
######### SUBSET TREES ###########
##################################

cd $dir/Simulations_folders/$path/$nameDir/
echo "Subset trees"
STARTTIME2=$(date +%s)

generations=$(seq 0 20 200)

cd $dir/Simulations_folders/$path/$nameDir/
for i in $(seq 1 1 $nbsimu)
do
	echo "python $dir/Subset_trees/subset_trees_villages_unilineal_descent.py -s $dir/Simulations_folders/$path/$nameDir/ -rep $i -g generations --sample-size $sample_size -K $K_total -og False -d $descent_rule -o $dir/Simulations_folders/$path/$nameDir/$i/ -t $dir/Tables/Simulations_metrics/$path/$nameDir > $i/outputPy${i}.txt"
done > launcher.txt
parallel -a launcher.txt -j $cores

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME2)) seconds to complete this task"

##################################
######## COMPUTE METRICS #########
##################################
echo "Compute metrics"
STARTTIME4=$(date +%s)
cd $dir/
rm -rf $dir/Tables/Pi/$path/$nameDir
mkdir -p $dir/Tables/Pi/$path/$nameDir
python Python_scripts/Pi_villages.py -s $dir/Simulations_folders/$path/$nameDir/ -ph $pseudohap -o $dir/Tables/Pi/$path/$nameDir/

python Python_scripts/Global_Pi_villages.py -s $dir/Simulations_folders/$path/$nameDir/ -ph $pseudohap -o $dir/Tables/Pi/$path/$nameDir/

cd $dir/
rm -rf $dir/Tables/$path/Fst/$nameDir/
mkdir -p $dir/Tables/Fst/$path/$nameDir/
python Python_scripts/Global_Fst_villages.py -s $dir/Simulations_folders/$path/$nameDir/ -d $descent -ph $pseudohap -o $dir/Tables/Fst/$path/$nameDir/

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME4)) seconds to complete this task"

##################################
########## WRITE NEXUS ###########
##################################
echo "Write .nex files"
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
	echo "python $dir/Subset_trees/write_nexus.py -s $dir/Simulations_folders/$path/$nameDir/ -rep $i --sample-size 20 -K $K_total -gen 200 -o $dir/BEAST/$path/$nameDir/$i/ > $dir/Simulations_folders/$path/$nameDir/$i/outputPyNex${i}.txt"
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
echo "Generate .xml files"
# create xml file
for i in $(seq 1 1 $nbsimu)
do
	echo "beastgen template_Mito.xml $path/$nameDir/$i/Sim_Mito.nex $path/$nameDir/$i/Sim_Mito.xml"
	echo "beastgen template_Y.xml $path/$nameDir/$i/Sim_Y.nex $path/$nameDir/$i/Sim_Y.xml"
done > launcher.txt
parallel -a launcher.txt -j $cores

cd $dir/BEAST/$path/$nameDir/
# beast inference
echo "Run BEAST"
for i in $(seq 1 1 $nbsimu)
do
	echo "beast -working -beagle_scaling always -overwrite $i/Sim_Mito.xml"
	echo "beast -working -beagle_scaling always -overwrite $i/Sim_Y.xml"
done > launcher.txt
parallel -a launcher.txt -j $cores

echo "Generate Skyline Plots"
for i in $(seq 1 1 $nbsimu)
do
	echo "Rscript $dir/Figures/R_scripts/skyline_plots.R '$dir/BEAST/$path/$nameDir/$i/'"
done > launcher.txt
parallel -a launcher.txt -j $cores

# Remove intermediate files
cd $dir/BEAST/$path/$nameDir/
for i in $(seq 1 1 $nbsimu)
do
	cd $i
	rm -v !(skyline*)
	cd $dir/BEAST/$path/$nameDir/
done

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME3)) seconds to complete this task"