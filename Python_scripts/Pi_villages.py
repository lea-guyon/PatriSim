import os
import numpy as np
from metrics import alleles_count, Pi, Theta
import glob
import argparse

def parse_args() :
	parser = argparse.ArgumentParser(description='Compute mean nucleotide diversity for each chromosome and each replicate')
	parser.add_argument('-s', '--source', dest='path_source', required=True, help='Source folder')
	parser.add_argument('-o', '--output', dest='output', help='Output folder path') 
	args = parser.parse_args()
	return args.path_source, args.output

path_source, output = parse_args()

os.chdir(output)
f_out1 = open("Pi_X_mean_by_rep.txt", 'w') # file containing means of Pi(X) per replicate
header=['mf/m', 'Gen', 'Pi', 'Theta', 'Nb_SNPs']
print('\t'.join(header), file=f_out1)

f_out2 = open("Pi_Y_mean_by_rep.txt", 'w') # file containing means of Pi(Y) per replicate
header=['mf/m', 'Gen', 'Pi', 'Theta', 'Nb_SNPs']
print('\t'.join(header), file=f_out2)

f_out3 = open("Pi_A_mean_by_rep.txt", 'w') # file containing means of Pi(A) per replicate
header=['mf/m', 'Gen', 'Pi', 'Theta', 'Nb_SNPs']
print('\t'.join(header), file=f_out3)

f_out4 = open("Pi_Mito_mean_by_rep.txt", 'w') # file containing means of Pi(Mito) per replicate
header=['mf/m', 'Gen', 'Pi', 'Theta', 'Nb_SNPs']
print('\t'.join(header), file=f_out4)

f_out5 = open("QPi_mean_by_rep.txt", 'w') # file containing means of QPi per replicate
header=['mf/m', 'Gen', 'Pi', 'Theta', 'Nb_SNPs']
print('\t'.join(header), file=f_out5)

f_out6 = open("Pi_mito_Y_mean_by_rep.txt", 'w') # file containing means of Pi(Y/Mito) per replicate
header=['mf/m', 'Gen', 'Pi', 'Theta', 'Nb_SNPs']
print('\t'.join(header), file=f_out6)

generations = [0, 20, 40, 60, 80, 100]

for rep in range (1, 201):
	for gen in generations :
		os.chdir(path_source)
		os.chdir(str(rep))

		# initialize indices
		pi_rep_X, theta_rep_X, nb_snps_rep_X = [], [], []
		pi_rep_Y, theta_rep_Y, nb_snps_rep_Y = [], [], []
		pi_rep_A, theta_rep_A, nb_snps_rep_A = [], [], []
		pi_rep_Mito, theta_rep_Mito, nb_snps_rep_Mito = [], [], []
		Q_pi_rep, Q_theta_rep, Q_nb_snps_rep = [], [], []
		pi_rep_mito_Y, theta_rep_mito_Y, nb_snps_rep_mito_Y = [], [], []

		chr_size_A = chr_size_X = chr_size_Y = 1e6
		chr_size_mito = 1e4

		filenames_A = sorted(glob.glob('Sim_A_*_{0}_gen_{1}_*.vcf'.format(rep, gen))) # list of files of replicate rep for autosomes
		if len(filenames_A) == 0 :
			continue
		filenames_X = sorted(glob.glob('Sim_X_*_{0}_gen_{1}_*.vcf'.format(rep, gen))) # list of files of replicate rep for X chr
		filenames_Y = sorted(glob.glob('Sim_Y_*_{0}_gen_{1}_*.vcf'.format(rep, gen))) # list of files of replicate rep for Y chr
		filenames_Mito = sorted(glob.glob('Sim_Mito_*_{0}_gen_{1}_*.vcf'.format(rep, gen))) # list of files of replicate rep for mt chr

		mf = filenames_A[0].split('_')[2]

		for j in range (0, len(filenames_A)):
			file_1A = filenames_A[j]
			AA1, AB1, nb_chr1A, nb_snps_A = alleles_count(file_1A) # Compute the number of alleles 0 and the number of alleles 1
			Pi_A = Pi(chr_size_A, AA1, AB1) # Compute Pi
			Theta_A = Theta(AA1, nb_chr1A, chr_size_A) # Compute Theta
			pi_rep_A += [Pi_A]
			theta_rep_A += [Theta_A]
			nb_snps_rep_A.append(nb_snps_A)

			file_1X = filenames_X[j]
			XA1, XB1, nb_chr1X, nb_snps_X = alleles_count(file_1X) # Compute the number of alleles 0 and the number of alleles 1
			Pi_X = Pi(chr_size_X, XA1, XB1) # Compute Pi
			Theta_X = Theta(XA1, nb_chr1X, chr_size_X) # Compute Theta
			pi_rep_X += [Pi_X]
			theta_rep_X += [Theta_X]
			nb_snps_rep_X.append(nb_snps_X)
			
			# Compute Q(pi)
			if Pi_A == 0 :
				Q_pi = float('nan')
			else :
				Q_pi = Pi_X/Pi_A

			# Compute Q(theta)
			if Theta_A == 0 :
				Q_theta = float('nan')
			else :
				Q_theta = Theta_X/Theta_A
			
			# Ratio of nb of SNPs
			if nb_snps_A == 0 :
				Q_nb_snps = float('nan')
			else :
				Q_nb_snps = nb_snps_X/nb_snps_A

			Q_pi_rep += [Q_pi]
			Q_theta_rep += [Q_theta]
			Q_nb_snps_rep.append(Q_nb_snps)

		for j in range (0, len(filenames_Y)):
			file_1Y = filenames_Y[j]
			YA1, YB1, nb_chr1Y, nb_snps_Y = alleles_count(file_1Y) # Compute the number of alleles 0 and the number of alleles 1
			Pi_Y = Pi(chr_size_Y, YA1, YB1) # Compute Pi
			Theta_Y = Theta(YA1, nb_chr1Y, chr_size_Y) # Compute theta
			pi_rep_Y += [Pi_Y]
			theta_rep_Y += [Theta_Y]
			nb_snps_rep_Y.append(nb_snps_Y)

		for j in range (0, len(filenames_Mito)):
			file_1M = filenames_Mito[j]
			MA1, MB1, nb_chr1M, nb_snps_mito = alleles_count(file_1M) # Compute the number of alleles 0 and the number of alleles 1
			Pi_M = Pi(chr_size_mito, MA1, MB1) # Compute Pi
			Theta_M = Theta(MA1, nb_chr1M, chr_size_mito) # Compute theta
			pi_rep_Mito += [Pi_M]
			theta_rep_Mito += [Theta_M]
			nb_snps_rep_Mito.append(nb_snps_mito)
	        
		pi_rep_mito_Y = np.nanmean(pi_rep_Mito)/np.nanmean(pi_rep_Y) # Compute Pi(Mito/Y)
		theta_rep_mito_Y = np.nanmean(theta_rep_Mito)/np.nanmean(theta_rep_Y) # Compute Theta(Mito/Y)
		nb_snps_rep_mito_Y = np.nanmean(nb_snps_rep_Mito)/np.nanmean(nb_snps_rep_Y) # Compute Tajima's D (Mito/Y)

		# Compute means
		mean_pi_A = np.nanmean(pi_rep_A)
		mean_theta_A = np.nanmean(theta_rep_A)
		mean_nb_snps_A = np.nanmean(nb_snps_rep_A)
		os.chdir(output)
		line = [str(mf), str(gen), str(mean_pi_A), str(mean_theta_A), str(mean_nb_snps_A)]
		print('\t'.join(line), file=f_out3)

		mean_pi_X = np.nanmean(pi_rep_X)
		mean_theta_X = np.nanmean(theta_rep_X)
		mean_nb_snps_X = np.nanmean(nb_snps_rep_X)
		os.chdir(output)
		line = [str(mf), str(gen), str(mean_pi_X), str(mean_theta_X), str(mean_nb_snps_X)]
		print('\t'.join(line), file=f_out1)

		mean_pi_Y = np.nanmean(pi_rep_Y)
		mean_theta_Y = np.nanmean(theta_rep_Y)
		mean_nb_snps_Y = np.nanmean(nb_snps_rep_Y)
		os.chdir(output)
		line = [str(mf), str(gen), str(mean_pi_Y), str(mean_theta_Y), str(mean_nb_snps_Y)]
		print('\t'.join(line), file=f_out2)

		mean_pi_Mito = np.nanmean(pi_rep_Mito)
		mean_theta_Mito = np.nanmean(theta_rep_Mito)
		mean_nb_snps_Mito = np.nanmean(nb_snps_rep_Mito)
		os.chdir(output)
		line = [str(mf), str(gen), str(mean_pi_Mito), str(mean_theta_Mito), str(mean_nb_snps_Mito)]
		print('\t'.join(line), file=f_out4)

		mean_pi_Q = np.nanmean(Q_pi_rep)
		mean_Q_theta = np.nanmean(Q_theta_rep)
		mean_Q_nb_snps = np.nanmean(Q_nb_snps_rep)
		os.chdir(output)
		line = [str(mf), str(gen), str(mean_pi_Q), str(mean_Q_theta), str(mean_Q_nb_snps)]
		print('\t'.join(line), file=f_out5)

		mean_pi_mito_Y = np.nanmean(pi_rep_mito_Y)
		mean_theta_mito_Y = np.nanmean(theta_rep_mito_Y)
		mean_nb_snps_mito_Y = np.nanmean(nb_snps_rep_mito_Y)
		os.chdir(output)
		line = [str(mf), str(gen), str(mean_pi_mito_Y), str(mean_theta_mito_Y), str(mean_nb_snps_mito_Y)]
		print('\t'.join(line), file=f_out6)

f_out1.close()
f_out2.close()
f_out3.close()
f_out4.close()
f_out5.close()
f_out6.close()