import os 
import numpy as np 
import pandas as pd
pd.plotting.register_matplotlib_converters()
import matplotlib.pyplot as plt 
import seaborn as sns 

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
#Run this if you want to see the entire sequence
#for sequence in SeqIO.parse('MN908947.fna','fasta'):
	#print(sequence.seq)
	#print(len(sequence),'Nucleotides')
#Load into a more easier-to-read-analyze format
#SeqIO.read() produces the basic information about the sequence
DNASeq = SeqIO.read('MN908947.fna','fasta')
print(DNASeq)
#The Coronavirus is an RNA sequence 
DNA = DNASeq.seq 
#Convert the DNA Sequence into mRNA sequence 
#Thymine replaced with Uracil
mRNA = DNA.transcribe()
#print(mRNA)
print('Size:',len(mRNA))

#Next, translate the RNA sequence into an Amino acid sequence 
#There is a translation table, which gives the NCBI genetic code number or name
#The '*' is a separator for proteins 
#Note that there are fewer sequences in the protein since 3 mRNAs are used to 
#produce a single subunit of protein
aminoAcid = mRNA.translate(table = 1,cds = False)
print('Amino Acid',aminoAcid)
print('Length of Protein:',len(aminoAcid))
print('Length of original RNA:',len(mRNA))
print(CodonTable.unambiguous_rna_by_name['Standard'])

#Using the standard Codon table, identify all the chains of amino acids (proteins)
#marked by *. Remove any sequence less than 20 amino acids long 
#This is the smallest known functional protein 
proteins = aminoAcid.split('*')
dfProteins = pd.DataFrame(proteins)
print(dfProteins.describe())
print(dfProteins[:10])
print("Total Proteins:",len(dfProteins))

def getLength(x):
	return len(x)

def getString(x):
	return ''.join(x)

dfProteins.rename(columns = {0:'Sequence'},inplace = True)
dfProteins['Length'] = dfProteins['Sequence'].apply(getLength)
funcProteins = dfProteins[dfProteins['Length']>=20]
funcProteins['SequenceString'] = funcProteins['Sequence'].apply(getString)
funcProteins.index = range(0,len(funcProteins))
funcProteins.info()
print("\n")
print(funcProteins.describe())
print(funcProteins[:10])

#Clear up some memory
del dfProteins

#Analyze the proteins 
from Bio.SeqUtils import ProtParam
poiList = []; mwList = []
for r in proteins[:]:
	X = ProtParam.ProteinAnalysis(str(r))
	POI = X.count_amino_acids()
	poiList.append(POI)
	molecularWeight = X.molecular_weight()
	mwList.append(molecularWeight)
	"""
	print("Protein of Interest:",POI)
	print("Percentage of Amino Acids:",
		X.get_amino_acids_percent())
	print("Molecular Weight:",
		mwList)
	print("Aromaticity:",
		X.aromaticity())
	print("Flexibility:",
		X.flexibility())
	print("Isoelectric point:",
		X.isoelectric_point())
	print("Secondary Structure Fraction:",
		X.secondary_structure_fraction())
	"""

#Plot the results 
mwDF = pd.DataFrame(data = mwList,columns = ['Molecular Weights'])
#Plot Proteins of interest 
poiDict = poiList[48]
plt.figure(figsize = (10,6))
plt.bar(poiDict.keys(),list(poiDict.values()),align = 'center')
plt.show()

#Compare this virus with SARS and MERS 
#Alignment using pairwise algorithm
from Bio import pairwise2
sars = SeqIO.read("","fasta")
mers = SeqIO.read("","fasta")
cov2 = SeqIO.read("","fasta")
sars_cov = pairwise2.align.globalxx(sars.seq, cov2.seq, one_alignment_only=True, score_only=True)
print('SARS/COV Similarity (%):', sars_cov / len(sars.seq) * 100)
sars_mers = pairwise2.align.globalxx(sars.seq, mers.seq, one_alignment_only=True, score_only=True)
print('SARS/MERS Similarity (%):', sars_mers / len(sars.seq) * 100)
cov2_mers = pairwise2.align.globalxx(cov2.seq, mers.seq, one_alignment_only=True, score_only=True)
print('MERS/COV Similarity (%):', mers_cov / len(mers.seq) * 100)









