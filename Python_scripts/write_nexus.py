import msprime, tskit
import os
import numpy as np
import random
import glob
import argparse

def parse_args() :
	parser = argparse.ArgumentParser(description='Split trees into A, X, Y and mito and generate vcf files')
	parser.add_argument('-s', '--source', dest='path_source', required=True, help='Source directory')
	parser.add_argument('-rep', '--replicat', dest='rep', type = int, required=True, help='Replicate number')
	parser.add_argument('--sample-size', dest='sample_size', type = int, required=True, help='Number (even) of individuals sampled per village')
	parser.add_argument('-K', '--carrying-capacity', dest='K', type = int, required=True, help='Total carrying capacity of the simulation')
	parser.add_argument('-gen', '--generation', dest='gen', type = int, required=True, help='Generation')
	parser.add_argument('-o', '--output', dest = 'output', required = True, help = 'Output file path')
	args = parser.parse_args()
	return args.path_source, args.rep, args.sample_size, args.K, args.gen, args.output

path_source, rep, sample_size, K, gen, output = parse_args()

###### Exceptions ######
if sample_size % 2 != 0 :
	raise Exception("Sample size isn't an even number")

chr_size = [1e6, 1e6, 1e4]

# change seed for each simulation
seed = str(random.sample(range(1,1000000000), 1)[0]) 

os.chdir(path_source)
os.chdir(str(rep))

###### Load file ######
filename = glob.glob('Sim_*{0}_gen_{1}.trees'.format(rep, gen))[0]
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

###### Sample individuals ######
indM, indF = [], []
M, F = [], []
nodes_M = []
nodes_F = []
villages = {} # associate village marker with individuals

for mut in ts.mutations() :
	if mut.site != 0 :
		continue
	if mut.time != 0 :
		continue
	mut_list = mut.metadata['mutation_list'] 
	marker = mut_list[0]['mutation_type']
	if marker not in villages :
		villages[marker] = []
	node = mut.node
	ind = ts.node(node).individual
	if ind not in villages[marker] :
		villages[marker].append(ind)

# sample N men and N women per village
N = int(sample_size/2)
for ind in ts.individuals() :
	if ind.metadata['sex'] == 1 :
		indM.append(ind.id)
	else :
		indF.append(ind.id)

for village in villages :
	vilM = [ind for ind in indM if ind in villages[village]]
	if len(vilM) < N :
		continue # not enough men
	M += random.sample(vilM, N)
	vilF = [ind for ind in indF if ind in villages[village]]
	if len(vilF) < N :
		continue # not enough women
	F += random.sample(vilF, N)

nodes, nodes_X, nodes_Y, nodes_mito = {}, {}, {}, {}
nodes_M, nodes_F = [], []
for ind in ts.individuals() :
	nodes_M.extend([node for node in ind.nodes if ind.id in M])
	nodes_F.extend([node for node in ind.nodes if ind.id in F])
nodes = nodes_M + nodes_F

nodes_ts_Y = [i for i in nodes_M if i in id_Y] # nodes' IDs in ts
nodes_Y = [ts_Y_map[1][i] for i in nodes_ts_Y] # nodes' IDs in recap_ts_Y

nodes_ts_X = [i for i in nodes if i in id_X] # nodes' IDs in ts
nodes_X = [ts_X_map[1][i] for i in nodes_ts_X] # nodes' IDs in recap_ts_X

nodes_ts_mito = [i for i in nodes_F if i in id_Mito] # nodes' IDs in ts
nodes_mito = [ts_mito_map[1][i] for i in nodes_ts_mito] # nodes' IDs in recap_ts_mito

###### Add mutations ######
mutated_ts_Y = msprime.sim_mutations(ts_Y, rate = 2.5e-8, random_seed = seed, model = "jc69", keep = False)
mutated_ts_mito = msprime.sim_mutations(ts_mito, rate = 5.5e-7, random_seed = seed, model = "jc69", keep = False)

###### Simplify ######
mutated_ts_Y = mutated_ts_Y.simplify(nodes_Y)
mutated_ts_mito = mutated_ts_mito.simplify(nodes_mito)

#reference_sequence_A = tskit.random_nucleotides(mutated_ts_A.sequence_length)
reference_sequence_Y = tskit.random_nucleotides(mutated_ts_Y.sequence_length)
reference_sequence_mito = tskit.random_nucleotides(mutated_ts_mito.sequence_length)

###### Output .nex files ######
os.chdir(output)
with open("Sim_Y_{0}_gen_{1}.nex".format(rep, gen), "w") as nexus_file:
	mutated_ts_Y.write_nexus(nexus_file, include_trees = False, reference_sequence = reference_sequence_Y)

with open("Sim_Mito_{0}_gen_{1}.nex".format(rep, gen), "w") as nexus_file:
	mutated_ts_mito.write_nexus(nexus_file, include_trees = False, reference_sequence = reference_sequence_mito)