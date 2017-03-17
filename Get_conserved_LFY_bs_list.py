
### ----------------------------
### LFY BS discovery program
### ----------------------------

'''
This program allows to identify LFY binding sites (BS).
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks)).
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import vcf
import numpy as np
from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type=str, default= "LFY_scores_matrix_19nucl")
parser.add_argument("--threshold", "-th",nargs='+',type=int,default= -23)
args = parser.parse_args()

#python get_interdistancesV2.py -fac "LFY_scores_matrix_19nucl" -th -23 

factorTranscription = args.factor
threshold = args.threshold
                    
###################Parameters we can change#################
	
if factorTranscription == "LFY_scores_matrix_19nucl" :
	FastaFile = "LFY_bound_sequences.fas" 
	MatrixFile = "LFY_scores_matrix_19nucl.txt" 
	dependencyFile = "interdependent_bases_matrix_for_LFY.txt"

####################### This is useful to calculate interdependances ###########

codigo = { 'ACC': 5,
           'ATG': 14, 'AAG': 2, 'AAA': 0, 'ATC': 13, 'AAC': 1, 'ATA': 12,
           'AGG': 10, 'CCT': 23, 'ACT': 7, 'AGC': 9, 'ACA': 4, 'AGA': 8,
           'CAT': 19, 'AAT': 3, 'ATT': 15, 'CTG': 30, 'CTA': 28,
           'CTC': 29, 'CAC': 17, 'ACG': 6,'CAA': 16, 'AGT': 11, 'CCA': 20,
           'CCG': 22, 'CCC': 21, 'TAT': 51, 'GGT': 43, 'TGT': 59, 'CGA': 24,
           'CAG': 18, 'CGC': 25, 'GAT': 35, 'CGG': 26, 'CTT': 31, 'TGC': 57,
           'GGG': 42, 'TAG': 50, 'GGA': 40, 'TAA': 48, 'GGC': 41, 'TAC': 49,
           'GAG': 34, 'TCG': 54, 'TTA': 60, 'GAC': 33, 'CGT': 27, 'TTT': 63,
           'TCA': 52, 'GCA': 36, 'GTA': 44, 'GCC': 37, 'GTC': 45, 'GCG': 38,
           'GTG': 46, 'TTC': 61, 'GTT': 47, 'GCT': 39, 'TGA': 56, 'TTG': 62,
           'TCC': 53, 'TGG': 58, 'GAA': 32, 'TCT': 55}
           
codigoi = { "A" : "T", "C" : "G", "G" : "C", "T" : "A"}

###############################################################################

def get_score_matrix(Mdata) :
	matScore = map(float,Mdata)
	lenMotif = len(matScore)/4
	return (matScore, lenMotif)

def get_dependency_matrix(dependencyFile) : 
	G = open(dependencyFile,"r")
	dependency_file_content = G.read().replace("\r","\n") + "\n"
	G.close()
	
	num2 = re.compile(r"(?<![\d.])[0-9]+(?![\d.])")
	position_dependency_matrix = num2.findall(dependency_file_content)
	position_dependency_matrix = map(int, position_dependency_matrix)
	
	dependency_matrix = num.findall(dependency_file_content)
	dependency_matrix = map(float, dependency_matrix)
	
	dependency_matrix_associated_with_position = []
	index1 = 0
	index2 = 3
	index3 = 0
	index4 = 64
	for i in range(0, 3):
		dependency_matrix_associated_with_position.append(position_dependency_matrix[index1:index2])
		dependency_matrix_associated_with_position.append(dependency_matrix[index3:index4])
		index1 = index1 + 3
		index2 = index2 + 3
		index3 = index3 + 64
		index4 = index4 + 64
		
	return(dependency_matrix_associated_with_position)
	
def divide(a, b):
    if b == 0:
        return np.nan
    else: 
        return a/b

def seq_c(site):
        site_i = site[-1::-1]
        site_c = ""
        for x in site_i:
		y = codigoi[x]
                site_c = site_c + y
        return site_c   

def add_scores_associated_with_interdependent_positions(dependency_matrix,scoreStrandPos,scoreStrandNeg,strandPos):
	cStrand = ""
	for lettre in seq_c(strandPos):
		cStrand = lettre + cStrand
	cStrand = cStrand[::-1]
	
	site1 = ""
	Csite1 = ""
	site2 = ""
	Csite2 = ""
	site3 = ""
	Csite3 = ""
	for i in dependency_matrix[0]:
		site1 = site1 + strandPos[i-1]
		Csite1 = Csite1 + cStrand[i-1]
	for i in dependency_matrix[2]:
		site2 = site2 + strandPos[i-1]
		Csite2 = Csite2 + cStrand[i-1]
	for i in dependency_matrix[4]:
		site3 = site3 + strandPos[i-1]
		Csite3 = Csite3 + cStrand[i-1]
	#print("dependency_matrix[1][0] : ",dependency_matrix[1][0])
	scoreStrandPos = scoreStrandPos + dependency_matrix[1][codigo[site1]] + dependency_matrix[3][codigo[site2]] + dependency_matrix[5][codigo[site3]]
	scoreStrandNeg = scoreStrandNeg + dependency_matrix[1][codigo[Csite1]] + dependency_matrix[3][codigo[Csite2]] + dependency_matrix[5][codigo[Csite3]]
	return(scoreStrandPos, scoreStrandNeg)

def get_list_LFY_binding_sites(matScore,matRev,FastaFile,dependency_matrix,threshold,factorTranscription):
	
	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

	print "  There are %s sequence(s) to analyze"%(len(sequences))
	
	list_of_the_LFY_binding_sites = []
	# We apply a loop on all the fasta sequences:
	for s in sequences:
		
		# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
		good_score_positions = [] 
		
		# This line allows to retrieve the DNA sequence
		seq = sequences[s].seq
		seq_id = sequences[s].id
		chrom = re.split(':',seq_id)
		pos = re.split(':|-',seq_id)
		
		# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
		for c in range(len(seq) - (lenMotif -1)):
			strandPos = seq[c:c+lenMotif].upper()
			test = 0
			for nu in strandPos :
				if nu not in ["A","C","G","T"]:
					test = 1
			if test == 1:
				score = "NA"
			else :
				#These lines allows to calculate a score for one sub-sequence
				index = 0
				scoreStrandPos = 0
				scoreStrandNeg = 0 
				while index < lenMotif:
					if strandPos[index] == 'A':
						scoreStrandPos = scoreStrandPos + matScore[index*4]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4]
					elif strandPos[index] == 'C':
						scoreStrandPos = scoreStrandPos + matScore[index*4+1]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4+1]
					elif strandPos[index] == 'G':
						scoreStrandPos = scoreStrandPos + matScore[index*4+2]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4+2]
					elif strandPos[index] == 'T':
						scoreStrandPos = scoreStrandPos + matScore[index*4+3]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4+3]			
					index += 1
				
				# This function allows to add scores that are associated with interdependent positions
				scoreStrandPos, scoreStrandNeg = add_scores_associated_with_interdependent_positions(dependency_matrix,scoreStrandPos,scoreStrandNeg,strandPos)
				
				#These lines allows to retrieve the chromosome and the positions where there is a predicted binding site (score above the threshold fixed by the user) . 
				if scoreStrandPos > threshold or scoreStrandNeg > threshold:
					list_of_the_LFY_binding_sites.append([chrom[0],int(pos[1]) + c + 1, int(pos[1]) + c + 1 + 19,str(strandPos[0:19])])

	return(list_of_the_LFY_binding_sites)

########################################### About the main matrix #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas...

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''
####################################################################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read().replace("\r","\n") + "\n"
F.close()

# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list
import re
num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)

matScore, lenMotif = get_score_matrix(Mdata)

'''The following line allows to produce the reversed matrix
if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
So we can calculate with this reverse matrix, the score of the complementary strand.
'''
matRev = list(reversed(matScore))

'''The LFY matrix gets 3 * 3 bases that are interdependent, so we have to retrieve in a list the interdependances:
dependency_matrix = [[4, 5, 6], [-2.482, -2.8048, -5.8493, -1.9992, -4.463, -5.8493, -6.0225, -6.0225, ...], 
			[9, 10, 11], [-4.3503, -5.1942, -4.3503, -5.1942, -5.1942, -5.1942, -3.434, -5.1942, -3.9448, ...], 
			[14, 15, 16], [-6.0225, -6.0225, -6.0225, -6.0225, -1.1954, -5.8493, -5.8493, -5.8493, -6.0225, -6.0225, -6.0225, ...]]
'''
dependency_matrix = get_dependency_matrix(dependencyFile)


########## get list[chromosome, first position of a BS, last position of a BS, 'sequence of the BS']:
list_of_the_LFY_binding_sites = get_list_LFY_binding_sites(matScore,matRev,FastaFile,dependency_matrix,threshold,factorTranscription)
print("list_of_the_LFY_binding_sites : ",list_of_the_LFY_binding_sites)

###### STEP 2: Applying a loop on all the vcf files and retrieve the new sequences

#dir = '/work/adrien/program_including_VCF_files/vcf_files'
#for root, dirs, filenames in os.walk(dir):
	#for f in filenames:
		#FILE = open(os.path.join(root, f), 'r')
		#vcf_reader = vcf.Reader(open('FILE', 'r'))
		#for record in vcf_reader:
			## check if the polymorphism is contained in the list_of_the_LFY_binding_sites
			## if this is the case, replace the original sequence by the new one.
			
# Redo the get_list_LFY_binding_sites function with the new list
# list_of_the_LFY_binding_sites = get_list_LFY_binding_sites(matScore,matRev,FastaFile,dependency_matrix,threshold,factorTranscription)

# You will only have the conserved LFY BS ! 

