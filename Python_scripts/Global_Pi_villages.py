import os
import numpy as np
from metrics import alleles_count_global, Pi
import re
import glob
import argparse

def parse_args() :
	parser = argparse.ArgumentParser(description='Script permettant des Pi globaux')
	parser.add_argument('-s', '--source', dest='path_source', required=True, help='Chemin vers le dossier source')
	parser.add_argument('-o', '--output', dest='output', help='Le chemin du fichier de sortie') 
	args = parser.parse_args()
	return args.path_source, args.output

path_source, output = parse_args()

os.chdir(output)
f_out = open("Global_Pi.txt", 'w') # file containing means of Pi(X) per replicate
header=['CHR', 'Gen', 'Pi']
print('\t'.join(header), file=f_out)

generations = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220]

for rep in range (1, 201):
	for gen in generations :
		os.chdir(path_source)
		os.chdir(str(rep))

		#initialize indices
		pi_rep_X = []
		pi_rep_Y = []
		pi_rep_A = []
		pi_rep_Mito = []
		Q_pi_rep = []
		pi_rep_mito_Y = []

		chr_size_A = chr_size_X = chr_size_Y = 1e6
		chr_size_mito = 1e4

		filenames_A = sorted(glob.glob('Sim_A_*{0}_gen_{1}_*[1-9]*.vcf'.format(rep, gen)))  # list of files of replicate rep for autosomes
		filenames_X = sorted(glob.glob('Sim_X_*{0}_gen_{1}_*[1-9]*.vcf'.format(rep, gen)))  # list of files of replicate rep for X chr
		filenames_Y = sorted(glob.glob('Sim_Y_*{0}_gen_{1}_*[1-9]*.vcf'.format(rep, gen)))  # list of files of replicate rep for Y chr
		filenames_Mito = sorted(glob.glob('Sim_Mito_*{0}_gen_{1}_*[1-9]*.vcf'.format(rep, gen)))  # list of files of replicate rep for mt chr

		if filenames_A == [] or filenames_X == [] or filenames_Y == [] or filenames_Mito == []:
			continue
		
		# autosomes
		A, B = alleles_count_global(filenames_A) # Compute the number of alleles 0 and the number of alleles 1
		Pi_A = Pi(chr_size_A, A, B)
		pi_rep_A += [Pi_A]

		os.chdir(output)
		line = ['A', str(gen), str(Pi_A)]
		print('\t'.join(line), file=f_out)

		# X chromosomes
		os.chdir(path_source)
		os.chdir(str(rep))
		A, B = alleles_count_global(filenames_X) # Compute the number of alleles 0 and the number of alleles 1
		Pi_X = Pi(chr_size_X, A, B) # Compute Pi
		pi_rep_X += [Pi_X]

		os.chdir(output)
		line = ['X', str(gen), str(Pi_X)]
		print('\t'.join(line), file=f_out)

		os.chdir(path_source)
		os.chdir(str(rep))

		# Y chromosomes
		A, B = alleles_count_global(filenames_Y) # Compute the number of alleles 0 and the number of alleles 1
		Pi_Y = Pi(chr_size_Y, A, B) # Compute Pi
		pi_rep_Y += [Pi_Y]

		os.chdir(output)
		line = ['Y', str(gen), str(Pi_Y)]
		print('\t'.join(line), file=f_out)

		# mt chromosomes
		os.chdir(path_source)
		os.chdir(str(rep))
		A, B = alleles_count_global(filenames_Mito) # Compute the number of alleles 0 and the number of alleles 1
		Pi_M = Pi(chr_size_mito, A, B) # Compute Pi
		pi_rep_Mito += [Pi_M]

		os.chdir(output)
		line = ['Mito', str(gen), str(Pi_M)]
		print('\t'.join(line), file=f_out)

		# Compute Q(pi)
		if Pi_A == 0 :
			Q_Pi = float('nan')
		else :
			Q_Pi = Pi_X/Pi_A
		Q_pi_rep += [Q_Pi]

		# Compute Pi(mt/Y)
		if Pi_Y == 0 :
			Pi_mito_Y = float('nan')
		else :
			Pi_mito_Y = Pi_M/Pi_Y
		pi_rep_mito_Y += [Pi_mito_Y]

	
		os.chdir(output)
		line = ['X/A', str(gen), str(Q_Pi)]
		print('\t'.join(line), file=f_out)

		os.chdir(output) 
		line = ['mito/Y', str(gen), str(Pi_mito_Y)]
		print('\t'.join(line), file=f_out)

f_out.close()
