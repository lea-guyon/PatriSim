import random

def alleles_count(vcf_file, pseudohap = False) :

	"""
    Description : Compute the number of alleles 0 and the number of alleles 1 
    for each SNP in a sampled population or village

	Arguments : vcf_file = VCF file for a given chromosome, a given replicate and a given generation
				chr_size = size of chromosomes
				pseudohap = boolean indicating if genomes are pseuhaploidized (default = False)

    Return : A = dictionnary associating each SNP with the number of alleles 0 in the sample
             B = dictionnary associating each SNP with the number of alleles 1 in the sample
             nb_chr = nb of chromosomes in the sample
	"""

	f1 = open(vcf_file, 'r')

	A = {}
	B = {}
	nb_chr = 0
	nb_snps = 0

	for ligne in f1:
        # skip header
		if ligne.startswith('#'):
			continue

		ligne = ligne[:-1].split('\t')
		line = ligne[9:] # select only SNP information

		if pseudohap == 'True' :
			coef_A = random.randint(0,1)
			coef_B = 1 - coef_A
			A[ligne[1]] = float(line.count('0|0') + coef_A * line.count('0|1') + coef_A * (line.count('1|0') + line.count('2|0') + line.count('0|2')) + line.count('0')) # count number of 0
			B[ligne[1]] = float(line.count('1|1') + line.count('2|2') + line.count('2|1') + line.count('1|2') + coef_B * (line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2')) + line.count('1') + line.count('2')) # count number of 1 (and 2)
			
			if nb_chr == 0 :
				nb_chr = A[ligne[1]] + B[ligne[1]]
				#nb_chr = float(2 * (line.count('0|0') + line.count('1|1') + line.count('2|2') + line.count('2|1') + line.count('1|2') + line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2')) + line.count('0') + line.count('1') + line.count('2'))
			
			if A[ligne[1]] != 0 and A[ligne[1]] != nb_chr :
				nb_snps += 1 # count number of polymorphic snps
		else :
			A[ligne[1]] = float(2 * line.count('0|0') + line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2') + line.count('0')) # count number of 0
			B[ligne[1]] = float(2 * line.count('1|1') + 2 * line.count('2|2') + 2 * line.count('2|1') + 2 * line.count('1|2') + line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2') + line.count('1') + line.count('2')) # count number of 1 (and 2)

			if nb_chr == 0 :
				nb_chr = A[ligne[1]] + B[ligne[1]] # nb of chromosomes

			if A[ligne[1]] != 0 and A[ligne[1]] != nb_chr :
				nb_snps += 1 # count number of polymorphic snps

	f1.close()
	return(A, B, nb_chr, nb_snps)

def alleles_count_global(filenames, pseudohap = False) :

	"""
    Description : Compute the number of alleles 0 and the number of alleles 1 
    for each SNP for all populations or villages

	Arguments : filenames = list of VCF files for a given chromosome, a given replicate and a given generation
				pseudohap = boolean indicating if genomes are pseuhaploidized (default = False)

	Return : 	A = dictionnary associating each file (ie sample) with a dictionnary associating each 
					SNP with the number of alleles 0 in the sample
				B = dictionnary associating each file (ie sample) with a dictionnary associating each 
					SNP with the number of alleles 1 in the sample
	"""

	A = {}
	B = {}
	for i in filenames :
		f = open(i, 'r')
		
		for ligne in f:
			# skip header
			if ligne[0] == '#':
				continue

			ligne = ligne[:-1].split('\t')
			line = ligne[9:]

			if ligne[1] not in A :
				A[ligne[1]] = 0
				B[ligne[1]] = 0

			if pseudohap == 'True' :
				coef_A = random.randint(0,1)
				coef_B = 1 - coef_A
				A[ligne[1]] += float(line.count('0|0') + coef_A * line.count('0|1') + coef_A * (line.count('1|0') + line.count('2|0') + line.count('0|2')) + line.count('0')) # count number of 0
				B[ligne[1]] += float(line.count('1|1') + line.count('2|2') + line.count('2|1') + line.count('1|2') + coef_B * (line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2')) + line.count('1') + line.count('2')) # count number of 1 (and 2)

			else :
				A[ligne[1]] += float(2 * line.count('0|0') + line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2') + line.count('0')) # count number of 0
				B[ligne[1]] += float(2 * line.count('1|1') + 2 * line.count('2|2') + 2 * line.count('2|1') + 2 * line.count('1|2') + line.count('0|1') + line.count('1|0') + line.count('2|0') + line.count('0|2') + line.count('1') + line.count('2')) # count number of 1 (and 2)

	f.close()
	return(A, B)

def Pi(chr_size, A, B) :

	"""
	Description : compute nucleotide diversity (ie expected heterozygosity)

	Arguments : chr_size = size (in nucleotides) of the modeled chromosomes
				A = dictionnary associating each SNP with the number of alleles 0 in the sample
				B = dictionnary associating each SNP with the number of alleles 1 in the sample
	
	Return : pi = mean nucleotide diversity over whole chromosomes
	"""
	Pi = 0
	for snp in A :
		if A[snp] + B[snp] <= 1 :
			continue
		Pi += A[snp]*B[snp]/(0.5*(A[snp]+B[snp]-1)*(A[snp]+B[snp]))
	pi = Pi/float(chr_size)

	return(pi)

def Theta(A, nb_chr, chr_size) :

	"""
	Description : compute Watterson's theta

	Arguments : A = dictionnary associating each SNP with the number of alleles 0 in the sample
				nb_chr = nb of chromosomes in the sample
				chr_size = size (in nucleotides) of the modeled chromosomes

	Return : theta = mean Watterson's theta over whole chromosomes
	"""
	A_copy = A.copy()
	for snp in A_copy :
		if A[snp] == 0 or A[snp] == nb_chr :
			del A[snp]
	S = len(A)  # nb of polymorphic sites
	n = int(nb_chr)
	a = sum([1/i for i in range(1, n)])
	if a == 0 :
		Theta = float('nan')
	else :
		Theta = S/a
	theta = Theta/chr_size
	return(theta)
