import msprime, pyslim, tskit
import os
import numpy as np
import random
import glob
import argparse
import warnings

def parse_args() :
	parser = argparse.ArgumentParser(description='Split trees into A, X, Y and mito and generate vcf files')
	parser.add_argument('-s', '--source', dest='path_source', required=True, help='Source directory')
	parser.add_argument('-rep', '--replicat', dest='rep', type = int, required=True, help='Replicate number')
	parser.add_argument('--sample-size', dest='sample_size', type = int, required=True, help='Number (even) of individuals sampled per village')
	parser.add_argument('-K', '--carrying-capacity', dest='K', type = int, required=True, help='Total carrying capacity of the simulation')
	parser.add_argument('-d', '--descent', dest='descent', required = True, help='Either "patrilineal" or "matrilineal"')
	parser.add_argument('-o', '--output', dest = 'output', required = True, help = 'Output file path')
	parser.add_argument('-t', '--output-table', dest = 'output_table', required = True, help = 'Output table path')
	args = parser.parse_args()
	return args.path_source, args.rep, args.sample_size, args.K, args.descent, args.output, args.output_table

path_source, rep, sample_size, K, descent, output, output_table = parse_args()

###### Exceptions ######
if sample_size % 2 != 0 :
	raise Exception("Sample size isn't an even number")

# ignore msprime warning for time units mismatch
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

generations = [0, 20, 40, 60, 80, 100, 220]
chr_size = [1e6, 1e6, 1e4]

# change seed for each simulation
seed = str(random.sample(range(1,1000000000), 1)[0]) 

for gen in generations : 
	os.chdir(path_source)
	os.chdir(str(rep))

	###### Load file ######
	print(glob.glob('Sim_*{0}_gen_{1}.trees'.format(rep, gen)))
	filename = glob.glob('Sim_*{0}_gen_{1}.trees'.format(rep, gen))[0]
	mf = filename.split('_')[1]
	ts = tskit.load(filename)

	###### Split trees into 3 for autosomes, X/Y and mtDNA ######
	ts_A = ts.keep_intervals(np.array([0, chr_size[0] + 1], ndmin=2))
	ts_A = ts_A.trim()  # remove empty intervals

	ts_X_Y = ts.keep_intervals(np.array([chr_size[0] + 1, chr_size[1] + chr_size[0] + 3], ndmin=2))
	ts_X_Y = ts_X_Y.trim()  # remove empty intervals

	ts_mito = ts.keep_intervals(np.array([chr_size[1] + chr_size[0] + 3, chr_size[2] + chr_size[1] + chr_size[0] + 5], ndmin=2))
	ts_mito = ts_mito.trim()  # remove empty intervals

	print("There are ", ts_X_Y.num_trees, " X/Y trees")

	###### Make lists of nodes for each chromosome type ######
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

	###### Recapitate if there is more than 1 root ######
	tsA_max_roots = max(t.num_roots for t in ts_A.trees())
	if tsA_max_roots > 1 :		
		if ts_A.num_populations > 99 :
			print("more than 99 subpopulations")
			demography = msprime.Demography()
			demography.add_population(name="p1", initial_size=K)
			for pop in ts_A.populations():
				if pop.id == 0 :
					continue
				name = pop.metadata['name']
				demography.add_population(name=name, initial_size=int(K/ts_A.num_populations))
				#demography.add_population(name=name, initial_size=int(K/(ts_A.num_populations * 2)))
				demography.add_mass_migration(time=gen, source=name, dest="p1", proportion=1)		
			ts_A = pyslim.recapitate(ts_A, recombination_rate = 1.1e-8, random_seed = seed, demography = demography)
		else :
			recap_ts_A = pyslim.recapitate(ts_A, recombination_rate = 1.1e-8, ancestral_Ne = K, random_seed = seed)

	tsX_max_roots = max(t.num_roots for t in ts_X.trees())
	if tsX_max_roots > 1 :
		if ts_X.num_populations > 99 :
			demography = msprime.Demography()
			demography.add_population(name="p1", initial_size=3/4*K)
			for pop in ts_X.populations():
				if pop.id == 0 :
					continue
				name = pop.metadata['name']
				demography.add_population(name=name, initial_size=int(3*K/(4*ts_X.num_populations)))
				#demography.add_population(name=name, initial_size=int(3*K/(8*ts_X.num_populations)))
				demography.add_mass_migration(time=gen, source=name, dest="p1", proportion=1)
			ts_X = pyslim.recapitate(ts_X, recombination_rate = 1e-8, random_seed = seed, demography = demography) 
		else :
			recap_ts_X = pyslim.recapitate(ts_X, recombination_rate = 1e-8, ancestral_Ne = 3/4*K, random_seed = seed)

	###### Sample individuals ######
	villages = {} # associate village marker with individuals
	nodes = []
	for tree in ts.trees() :
		for mut in tree.mutations() :
			if mut.site != 0 :
				continue
			if mut.time != 0 :
				continue
			mut_list = mut.metadata['mutation_list'] 
			marker = mut_list[0]['mutation_type']
			nodes = [i for i in tree.leaves(mut.node)]
			if marker not in villages :
				villages[marker] = []
			villages[marker].extend(nodes)

	for marker in villages :
		villages[marker] = list(set(villages[marker]))

	nodes, nodes_X, nodes_Y, nodes_mito = {}, {}, {}, {}
	M, F = {}, {} # dictionnary of male (M) and female (F) samples for each village
	# sample N men and N women per village
	N = int(sample_size/2)
	indM, indF = [], []
	for ind in ts.individuals() :
		if ind.metadata['sex'] == 1 :
			indM.append(ind)
		else :
			indF.append(ind)

	for village in villages :
		print(village)
		vilM = [ind.id for ind in indM if ind.nodes[0] in villages[village]]
		print(len(vilM))
		if len(vilM) < N :
			continue # not enough men
		M[village] = random.sample(vilM, N)
		vilF = [ind.id for ind in indF if ind.nodes[0] in villages[village]]
		if len(vilF) < N :
			continue # not enough women
		F[village] = random.sample(vilF, N)

	for village in villages :
		nodes_M, nodes_F = [], []
		for ind in ts.individuals() :
			nodes_M.extend([node for node in ind.nodes if ind.id in M[village]])
			nodes_F.extend([node for node in ind.nodes if ind.id in F[village]])
		nodes[village] = nodes_M + nodes_F

		nodes_ts_Y = [i for i in nodes_M if i in id_Y] # nodes' IDs in ts
		nodes_Y[village] = [ts_Y_map[1][i] for i in nodes_ts_Y] # nodes' IDs in recap_ts_Y

		nodes_ts_X = [i for i in nodes[village] if i in id_X] # nodes' IDs in ts
		nodes_X[village] = [ts_X_map[1][i] for i in nodes_ts_X] # nodes' IDs in recap_ts_X

		nodes_ts_mito = [i for i in nodes_F if i in id_Mito] # nodes' IDs in ts
		nodes_mito[village] = [ts_mito_map[1][i] for i in nodes_ts_mito] # nodes' IDs in recap_ts_mito
    
	village_name = {4 : 1, 5 : 2, 6 : 3, 7 : 4, 8 : 5}

	###### Add mutations ######
	model = msprime.SLiMMutationModel(type = 1)

	mutated_ts_A = msprime.sim_mutations(ts_A, rate = 2.5e-8, random_seed = seed, model = model, keep = False)
	mutated_ts_X = msprime.sim_mutations(ts_X, rate = 2e-8, random_seed = seed, model = model, keep = False)
	mutated_ts_Y = msprime.sim_mutations(ts_Y, rate = 2.5e-8, random_seed = seed, model = model, keep = False)
	mutated_ts_mito = msprime.sim_mutations(ts_mito, rate = 5.5e-7, random_seed = seed, model = model, keep = False)
	
	###### Output VCF files ######
	for village in village_name :
		# select individuals belonging to the village
		ind_A = M[village] + F[village]
		with open("Sim_A_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village_name[village]), "w") as vcf_file:
			mutated_ts_A.write_vcf(vcf_file, contig_id = 'A', individuals = ind_A)

		ind_X = M[village] + F[village]
		with open("Sim_X_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village_name[village]), "w") as vcf_file:
			mutated_ts_X.write_vcf(vcf_file, contig_id = 'X', individuals = ind_X)

		ind_Y = []
		for indiv in mutated_ts_Y.individuals() :
			if list(indiv.nodes)[0] in nodes_Y[village] :
				ind_Y.append(indiv.id)
		with open("Sim_Y_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village_name[village]), "w") as vcf_file:
			mutated_ts_Y.write_vcf(vcf_file, contig_id = 'Y', individuals = ind_Y)

		ind_mito = []
		for indiv in mutated_ts_mito.individuals() :
			if list(indiv.nodes)[0] in nodes_mito[village] :
				ind_mito.append(indiv.id)
		with open("Sim_Mito_{0}_{1}_gen_{2}_village_{3}.vcf".format(mf, rep, gen, village_name[village]), "w") as vcf_file:
			mutated_ts_mito.write_vcf(vcf_file, contig_id = 'Mito', individuals = ind_mito)