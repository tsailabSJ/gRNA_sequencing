import sys
import pandas as pd
from Bio import SeqIO

def read_fasta_to_list(f):
	myList = []
	for r in SeqIO.parse(f, "fasta"):
		myList.append(str(r.seq).upper())
	return myList	
def revcomp(seq):
	try: ## python2
		tab = string.maketrans(b"ACTG", b"TGAC")
	except:  ## python3
		tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]
input=sys.argv[1]
gRNA_list = sys.argv[2]



def is_gRNA(seq,gRNA):
	if gRNA in seq:
		return 1
	if revcomp(gRNA) in seq:
		return 1
	return 0

def get_value(input,gRNA_list):
	df = pd.DataFrame(read_fasta_to_list(input))
	gRNA = pd.read_csv(gRNA_list,sep="\t",header=None)
	for i,r in gRNA.iterrows():
		df[r[0]] = [is_gRNA(seq,r[1]) for seq in df[0]]

	total = df.shape[0]
	df2 = df.copy()
	df2 = df2.drop([0],axis=1)
	df3 = pd.DataFrame((df2.sum()/total)*100)
	df3.columns = [f'%Read_{input}']
	return df3
