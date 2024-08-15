import sys
import os
import argparse
import gzip
import itertools
import subprocess
import gzip
import pandas as pd
import glob
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide as pa
from joblib import Parallel, delayed
from Levenshtein import distance
from Bio import pairwise2

def read_fastq_to_dict(f):
	## added collapse same sequences
	out = []
	# f1 = gzip.open(f, 'rt')
	f1 = open(f, 'rt')
	line1 = f1.readline()
	while (line1):
		first_line1 = line1.strip()
		# second_line1 = f1.readline().strip().replace("N","A")
		second_line1 = f1.readline().strip() # we will do replacement later
		third_line1 = f1.readline().strip()
		forth_line1 = f1.readline().strip()
		out.append([first_line1.split()[0],second_line1,third_line1,forth_line1])
		line1 = f1.readline()
	df = pd.DataFrame(out)
	df[4] = df.groupby([1])[0].transform('count')
	df = df.drop_duplicates(1)
	df.index = df[0]+"_"+df[4].astype(str)
	# df = df.head(n=100)
	return df[1].to_dict()
def revcomp(seq):
	try: ## python2
		tab = string.maketrans(b"ACTG", b"TGAC")
	except:  ## python3
		tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]



def extract_aligned_gRNA2(true_gRNA_sequence,s1,s2,all_N_index): # s1 is gRNA, s2 is read
	"""return gRNA sequence alignment in s1 and s2
	
	all_N_index to correct all N
	"""
	s1=str(s1)
	s2=str(s2)
	mutated_index = []
	if len(all_N_index)>0:
		count = -1
		for i in range(len(s2)):
			if s2[i]!="-":
				count += 1
				if count in all_N_index:
					mutated_index.append(i)
		s2 = list(s2)
		for i in mutated_index:
			s2[i]="N"
		s2 = "".join(s2)
	first_base = true_gRNA_sequence[0]
	last_base = true_gRNA_sequence[-1]
	first_base_pos = s1.index(first_base)
	last_base_pos = s1.rfind(last_base)+1
	return first_base_pos,s1[first_base_pos:last_base_pos],s2[first_base_pos:last_base_pos]

# true_gRNA_sequence = "CACCTACCTAAGAACCATCC"
# s1="---CACCTACC-TAAGAACCATCC---"
# s2="GCGCACCTACCCTAAGAACCATCCGTT"
# extract_aligned_gRNA2(true_gRNA_sequence,s1,s2)
def find_N(s,ch="N"):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def align_gRNA_and_scaffold(read,gRNA,scaffold):
	"""align gRNA sequences given the read
	
	
	output
	------
	
	user_gRNA: just the input gRNA sequence
	
	actual_gRNA: aligned gRNA sequence
	
	gRNA_distance: distance between user_gRNA and actual_gRNA
	
	stats
	
	added aligned gRNA sequence and scaffold sequence
	
	
	"""
	# scaffold="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
	is_perfect_gRNA = gRNA in read
	gRNA_distance = -1
	scaffold_distance = -1
	gRNA_start = -1
	obs_gRNA = "" # aligned gRNA sequence in read
	obs_scaffold = "" # aligned scaffold sequence in read
	aligned_gRNA = "" # aligned gRNA sequence of the given gRNA sequence, most likely to be the same, may have -
	aligned_scaffold = "" # aligned scaffold sequence of the given scaffold sequence
	status=""
	all_N_index = []
	# forward read
	#---- finding gRNA -----------------
	if is_perfect_gRNA:
		gRNA_start = read.index(gRNA)
		obs_gRNA = gRNA
		aligned_gRNA = gRNA
		gRNA_distance=0
	else:
		if "N" in read:
			status="containedN"
			# read = read.replace("N","A")
			all_N_index = find_N(read)
		s1=DNA(gRNA)
		s2=DNA(read.replace("N","A"))
		alignment, score, start_end_positions = pa(s1,s2)
		if score>=10:
			gRNA_start,aligned_gRNA,obs_gRNA = extract_aligned_gRNA2(gRNA,alignment[0],alignment[1],all_N_index)
			gRNA_distance = distance(gRNA,obs_gRNA.replace("-",""))
		else: # Failed to find gRNA, try to find scaffold and then infer gRNA
			s1=DNA(scaffold)
			alignment, score, start_end_positions = pa(s1,s2,match_score=2,gap_open_penalty=4)
			if score>=10:
				scaffold_start,aligned_scaffold,obs_scaffold = extract_aligned_gRNA2(scaffold,alignment[0],alignment[1],all_N_index)
				gRNA_end = scaffold_start
				gRNA_start = max(3,gRNA_end - len(gRNA))
				obs_gRNA = read[gRNA_start:gRNA_end]
				gRNA_distance = distance(gRNA,obs_gRNA)
				obs_scaffold = read[scaffold_start:scaffold_start+len(scaffold)]
				scaffold_distance = distance(scaffold, obs_scaffold)
				alignments = pairwise2.align.globalms(scaffold, obs_scaffold,match=2, mismatch=-1, open=-5,extend=-2)
				# pairwise2.align.globalxx(scaffold, obs_scaffold,match=2, mismatch=-1, open=-5,extend=-2)
				aligned_scaffold,obs_scaffold =alignments[0].seqA,alignments[0].seqB
				return [gRNA_start,gRNA,obs_gRNA,aligned_gRNA,gRNA_distance,scaffold,obs_scaffold,aligned_scaffold,scaffold_distance,status+"_+_"+"NoGrna"]
	#----------- finding scaffold -----------
	# gRNA is found, let's align scaffold
	if gRNA_start!=-1: 
		obs_gRNA_length = len(obs_gRNA.replace("-",""))
		obs_scaffold = read[gRNA_start+obs_gRNA_length:gRNA_start+obs_gRNA_length+len(scaffold)]
		if obs_scaffold == "":
			obs_scaffold = "N"
		scaffold_distance = distance(scaffold, obs_scaffold)
		alignments = pairwise2.align.globalms(scaffold, obs_scaffold,match=2, mismatch=-1, open=-5,extend=-2)
		try:
			aligned_scaffold,obs_scaffold =alignments[0].seqA,alignments[0].seqB
		except:
			print ("scaffold",scaffold)
			print ("obs_scaffold",obs_scaffold)
			print ("read",read)
			print (alignments)
		return [gRNA_start,gRNA,obs_gRNA,aligned_gRNA,gRNA_distance,scaffold,obs_scaffold,aligned_scaffold,scaffold_distance,"+"]
	# check reverse strand
	read_rv=revcomp(read)
	is_perfect_gRNA = gRNA in read_rv
	if is_perfect_gRNA:
		gRNA_start = read_rv.index(gRNA)
		obs_gRNA = gRNA
		aligned_gRNA = gRNA
		gRNA_distance=0
	else:
		if "N" in read_rv:
			all_N_index = find_N(read_rv)
		s1=DNA(gRNA)
		s2=DNA(read_rv.replace("N","A"))
		alignment, score, start_end_positions = pa(s1,s2)
		if score>10:
			# gRNA_start,obs_gRNA,gRNA_distance = extract_aligned_gRNA(alignment[0],alignment[1])
			gRNA_start,aligned_gRNA,obs_gRNA = extract_aligned_gRNA2(gRNA,alignment[0],alignment[1],all_N_index)
			gRNA_distance = distance(gRNA,obs_gRNA)
		else: # find scaffold
			s1=DNA(scaffold)
			alignment, score, start_end_positions = pa(s1,s2,match_score=2,gap_open_penalty=4)
			if score>=10:
				scaffold_start,aligned_scaffold,obs_scaffold = extract_aligned_gRNA2(scaffold,alignment[0],alignment[1],all_N_index)
				gRNA_end = scaffold_start
				gRNA_start = max(3,gRNA_end - len(gRNA))
				obs_gRNA = read_rv[gRNA_start:gRNA_end]
				gRNA_distance = distance(gRNA,obs_gRNA)
				obs_scaffold = read_rv[scaffold_start:scaffold_start+len(scaffold)]
				scaffold_distance = distance(scaffold, obs_scaffold)
				alignments = pairwise2.align.globalms(scaffold, obs_scaffold,match=2, mismatch=-1, open=-5,extend=-2)
				aligned_scaffold,obs_scaffold =alignments[0].seqA,alignments[0].seqB
				return [gRNA_start,gRNA,obs_gRNA,aligned_gRNA,gRNA_distance,scaffold,obs_scaffold,aligned_scaffold,scaffold_distance,status+"_-_"+"NoGrna"]
	if gRNA_start!=-1:
		obs_gRNA_length = len(obs_gRNA.replace("-",""))
		obs_scaffold = read_rv[gRNA_start+obs_gRNA_length:gRNA_start+obs_gRNA_length+len(scaffold)]
		if obs_scaffold == "":
			obs_scaffold = "N"
		scaffold_distance = distance(obs_scaffold,scaffold)
		alignments = pairwise2.align.globalms(scaffold, obs_scaffold,match=2, mismatch=-1, open=-5,extend=-2)
		aligned_scaffold,obs_scaffold =alignments[0].seqA,alignments[0].seqB
		return [gRNA_start,gRNA,obs_gRNA,aligned_gRNA,gRNA_distance,scaffold,obs_scaffold,aligned_scaffold,scaffold_distance,"-"]
	# failed to find gRNA, need to infer gRNA sequence
	gRNA_start=3
	obs_gRNA = read[gRNA_start:gRNA_start+len(gRNA)]
	obs_scaffold = read[gRNA_start+len(gRNA):gRNA_start+len(gRNA)+len(scaffold)]
	gRNA_distance = distance(obs_gRNA,gRNA)
	scaffold_distance = distance(obs_scaffold,scaffold)
	return [gRNA_start,gRNA,obs_gRNA,aligned_gRNA,gRNA_distance,scaffold,obs_scaffold,aligned_scaffold,scaffold_distance,"failed"]


file = sys.argv[1]
gRNA=sys.argv[2]
scaffold=sys.argv[3]
myDict = read_fastq_to_dict(file)
read_list = [*myDict]
print ("total reads",len(read_list))
values = [myDict[r] for r in read_list]
# out = Parallel(n_jobs=10,verbose=1,backend="loky")(delayed(is_scaffold_is_gRNA)(i,gRNA) for i in values)
out = Parallel(n_jobs=20,verbose=10,backend="multiprocessing")(delayed(align_gRNA_and_scaffold)(i,gRNA,scaffold) for i in values)
df = pd.DataFrame(out)
df.index = read_list
df.columns = ["gRNA_start","gRNA","obs_gRNA",'aligned_gRNA',"gRNA_distance",'scaffold',"obs_scaffold",'aligned_scaffold',"scaffold_distance","status"]
df['seq'] = values
df.to_csv(file.split(".fastq")[0]+".alignment.csv")










