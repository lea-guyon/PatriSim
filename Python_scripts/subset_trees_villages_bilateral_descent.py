import msprime, pyslim, tskit
import os
import numpy as np
import random
import glob
import warnings
import argparse

def parse_args() :
	parser = argparse.ArgumentParser(description='Split trees into A, X, Y and mito and generate vcf files')
	parser.add_argument('-s', '--source', dest='path_source', required=True, help='Source directory')
	parser.add_argument('-rep', '--replicat', dest='rep', type = int, required=True, help='Replicate number')
	parser.add_argument('-g', dest='generations', type=list, resuired=True, help='List of generation times when VCF files are output')
	parser.add_argument('--sample-size', dest='sample_size', type = int, required=True, help='Number (even) of individuals sampled per village')
	parser.add_argument('-K', '--carrying-capacity', dest='K', type = int, required=True, help='Total carrying capacity of the simulation')
	parser.add_argument('-o', '--output', dest = 'output', required = True, help = 'Output file path')
	args = parser.parse_args()
	return args.path_source, args.rep, args.g, args.sample_size, args.K, args.output

path_source, rep, generations, sample_size, K, output = parse_args()

# Exceptions
if sample_size % 2 != 0 :
	raise Exception("Sample size isn't an even number")

# ignore msprime warning for time units mismatch
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

for gen in generations : 
	os.chdir(path_source)
	os.chdir(str(rep))
	filename = glob.glob('Sim_*{0}_gen_{1}.trees'.format(rep, gen))[0]
	mf = filename.split('_')[1]
	ts = tskit.load(filename)

	# split trees into 3 for autosomes, X/Y and mito
	ts_A = ts.keep_intervals(np.array([0, 1e6+1], ndmin=2))
	ts_A = ts_A.trim()  # remove empty intervals

	ts_X_Y = ts.keep_intervals(np.array([1e6+1, 2e6+3], ndmin=2))
	ts_X_Y = ts_X_Y.trim()  # remove empty intervals

	ts_mito = ts.keep_intervals(np.array([2e6+3, 2e6+1e4+5], ndmin=2))
	ts_mito = ts_mito.trim()  # remove empty intervals

	print("There are ", ts_X_Y.num_trees, " X/Y trees")

	id_X_Y = [node.id for node in ts_X_Y.nodes() if node.time == 0]
	id_Y = []
	id_Mito = []

	for tree in ts_X_Y.trees():
		for mut in tree.mutations():
			id_Y += [i for i in tree.leaves(mut.node)] # list of nodes' ids for the Y chr
	id_Y = list(set(id_Y))
	id_X = [j for j in id_X_Y if j not in id_Y] # list of nodes' ids for the X chr

	for tree in ts_mito.trees():
		for mut in tree.mutations():
			id_Mito += [i for i in tree.leaves(mut.node)] # list of nodes' ids carrying a mutation representing mt chr
	id_Mito = list(set(id_Mito))

	ts_Y_map = ts_X_Y.simplify(id_Y, map_nodes=True, keep_input_roots=True) # Y chr tree + correspondances with the ids of ts
	ts_Y = ts_Y_map[0]

	ts_X_map = ts_X_Y.simplify(id_X, map_nodes=True, keep_input_roots=True) # X chr tree + correspondances with the ids of ts
	ts_X = ts_X_map[0]

	ts_mito_map = ts_mito.simplify(id_Mito, map_nodes=True, keep_input_roots=True) # mt chr tree + correspondances with the ids of ts
	ts_mito = ts_mito_map[0]

	print("There are ", ts_Y.num_trees, " Y trees")
	#for tree in ts_Y.trees():
	#	print(tree.draw(format="unicode"))
		
	print("There are", ts_X.num_trees, " X trees")
	#for tree in ts_X.trees():
	#	print(tree.draw(format="unicode"))

	print("There are", ts_mito.num_trees, " Mito trees")
	#for tree in ts_mito.trees():
	#	print(tree.draw(format="unicode"))

	seed = str(random.sample(range(1,1000000000), 1)[0]) # change seed for each simulation

	###### Recapitate if there is more than 1 root ######
	tsA_max_roots = max(t.num_roots for t in ts_A.trees())
	if tsA_max_roots > 1 :
		demography = msprime.Demography()
		demography.add_population(name="p1", initial_size=K)
		for pop in ts_A.populations():
			if pop.id == 0 :
				continue
			name = pop.metadata['name']
			demography.add_population(name=name, initial_size=int(K/(ts_A.num_populations)))
			demography.add_mass_migration(time=gen, source=name, dest="p1", proportion=1)
		ts_A = pyslim.recapitate(ts_A, recombination_rate = 1.1e-8, random_seed = seed, demography = demography)  

	tsX_max_roots = max(t.num_roots for t in ts_X.trees())
	if tsX_max_roots > 1 :
		demography = msprime.Demography()
		demography.add_population(name="p1", initial_size=3/4*K)
		for pop in ts_X.populations():
			if pop.id == 0 :
				continue
			name = pop.metadata['name']
			demography.add_population(name=name, initial_size=int(3*K/(4*ts_X.num_populations)))
			demography.add_mass_migration(time=gen, source=name, dest="p1", proportion=1)
		ts_X = pyslim.recapitate(ts_X, recombination_rate = 1e-8, random_seed = seed, demography = demography) 

	# sample individuals and simplify
	nodes_M = []
	nodes_F = []

	villages = set(list([ind.metadata['subpopulation'] for ind in ts_A.individuals()]))
	for village in villages :
		# select individuals belonging to the subpopulation 
		indF, indM = [], [] # initialize list of women and list of men
		for ind in ts.individuals() :
			if ind.metadata['subpopulation'] == village :
				if ind.metadata['sex'] == 0 :
					indF += [ind.id]
				else :
					indM += [ind.id]
		print(village, len(indM))
		# sample N men and N women per village
		N = int(sample_size/2)
		if len(indM) < N or len(indF) < N :
			continue # not enough women or not enough men
		M = random.sample(indM, N) # sample men
		F = random.sample(indF, N) # sample women
		for ind in ts.individuals() :
			nodes_M += [node for node in ind.nodes if ind.id in M]
			nodes_F += [node for node in ind.nodes if ind.id in F]
	nodes = nodes_M + nodes_F
	
	nodes_ts_Y = [i for i in nodes_M if i in id_Y] # nodes' IDs in ts
	nodes_Y = [ts_Y_map[1][i] for i in nodes_ts_Y] # nodes' IDs in recap_ts_Y

	nodes_ts_X = [i for i in nodes if i in id_X] # nodes' IDs in ts
	nodes_X = [ts_X_map[1][i] for i in nodes_ts_X] # nodes' IDs in recap_ts_X

	nodes_ts_mito = [i for i in nodes_F if i in id_Mito] # nodes' IDs in ts
	nodes_mito = [ts_mito_map[1][i] for i in nodes_ts_mito] # nodes' IDs in recap_ts_mito

	ind_A, ind_X, ind_Y, ind_mito = [], [], [], []

	for ind in ts_A.individuals() :
		if list(ind.nodes)[0] in nodes :
			ind_A.append(ind)

	for ind in ts_X.individuals() :
		if list(ind.nodes)[0] in nodes_X :
			ind_X.append(ind)

	for ind in ts_Y.individuals() :
		if list(ind.nodes)[0] in nodes_Y :
			ind_Y.append(ind)

	for ind in ts_mito.individuals() :
		if list(ind.nodes)[0] in nodes_mito :
			ind_mito.append(ind)

	# Add mutations
	model = msprime.SLiMMutationModel(type = 1)

	mutated_ts_A = msprime.sim_mutations(ts_A, rate = 2.5e-8, random_seed = seed, model = model, keep = False)
	mutated_ts_X = msprime.sim_mutations(ts_X, rate = 2e-8, random_seed = seed, model = model, keep = False)
	mutated_ts_Y = msprime.sim_mutations(ts_Y, rate = 2.5e-8, random_seed = seed, model = model, keep = False)
	mutated_ts_mito = msprime.sim_mutations(ts_mito, rate = 5.5e-7, random_seed = seed, model = model, keep = False)

	# Output VCFs 
	# for each subpopulation
	for village in villages :
		# select individuals belonging to the subpopulation (use _recap file because using mslim alters 'metadata')
		indiv_A = [ind.id for ind in ind_A if ind.metadata['subpopulation'] == village]
		with open(str(output) + "Sim_A_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village), "w") as vcf_file:
			mutated_ts_A.write_vcf(vcf_file, contig_id = 'A', individuals = indiv_A)

		# idem for X, Y and mito
		indiv_X = [ind.id for ind in ind_X if ind.metadata['subpopulation'] == village]
		with open(str(output) + "Sim_X_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village), "w") as vcf_file:
			mutated_ts_X.write_vcf(vcf_file, contig_id = 'X', individuals = indiv_X)

		indiv_Y = [ind.id for ind in ind_Y if ind.metadata['subpopulation'] == village]
		with open(str(output) + "Sim_Y_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village), "w") as vcf_file:
			mutated_ts_Y.write_vcf(vcf_file, contig_id = 'Y', individuals = indiv_Y)

		indiv_mito = [ind.id for ind in ind_mito if ind.metadata['subpopulation'] == village]
		with open(str(output) + "Sim_Mito_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village), "w") as vcf_file:
			mutated_ts_mito.write_vcf(vcf_file, contig_id = 'Mito', individuals = indiv_mito)

