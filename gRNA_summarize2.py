import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
import logomaker
from Bio import motifs
from Bio.Seq import Seq
import numpy as np
from joblib import Parallel, delayed
import sys
import os, io, random
import string
import numpy as np
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide as pa
from joblib import Parallel, delayed
from Levenshtein import distance
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from matplotlib.patches import Rectangle

# import panel as pn
# import panel.widgets as pnw
# pn.extension()

from bokeh.plotting import figure,save
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import export_png
from itertools import tee, count
from collections import Counter # Counter counts the number of occurrences of each item

def uniquify(seq, suffs = ('.0%s'%(x) for x in range(1, 1000))):
    """Make all the items unique by adding a suffix (1, 2, etc).

    `seq` is mutable sequence of strings.
    `suffs` is an optional alternative suffix iterable.
    """
    not_unique = [k for k,v in Counter(seq).items() if v>1] # so we have: ['name', 'zip']
    # suffix generator dict - e.g., {'name': <my_gen>, 'zip': <my_gen>}
    suff_gens = dict(zip(not_unique, tee(suffs, len(not_unique))))  
    for idx,s in enumerate(seq):
        try:
            suffix = str(next(suff_gens[s]))
        except KeyError:
            # s was unique
            continue
        else:
            # print (idx, suffix)
            seq[idx] += float(suffix)

def muscle_alignment(input,output):
	os.system("muscle -align %s -output %s"%(input,output))
	align = AlignIO.read(output, 'fasta')
	#os.system("rm %s"%(output))
	return align
def view_alignment(aln, fontsize="9pt", plot_width=800):
	"""Bokeh sequence alignment view"""

	#make sequence and id lists from the aln object
	aln.sort(key=lambda x:int(float(x.id)),reverse=False)
	seqs = [rec.seq for rec in (aln)]
	ids = [str(float(rec.id)) for rec in aln]	
	text = [i for s in list(seqs) for i in s]
	colors = get_colors(seqs)	
	N = len(seqs[0])
	S = len(seqs)	
	width = .4

	x = np.arange(1,N+1)
	y = np.arange(0,S,1)
	#creates a 2D grid of coords from the 1D arrays
	xx, yy = np.meshgrid(x, y)
	#flattens the arrays
	gx = xx.ravel()
	gy = yy.flatten()
	#use recty for rect coords with an offset
	recty = gy+.5
	h= 1/S
	#now we can create the ColumnDataSource with all the arrays
	source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
	plot_height = len(seqs)*15+50
	x_range = Range1d(0,N+1, bounds='auto')
	if N>100:
		viewlen=100
	else:
		viewlen=N
	#view_range is for the close up view
	view_range = (0,viewlen+1)
	tools="xpan, xwheel_zoom, reset, save"

	#entire sequence view (no text, with zoom)
	p = figure(title=None, plot_width= plot_width, plot_height=50,
			   x_range=x_range, y_range=(0,S), tools=tools,
			   min_border=0, toolbar_location='below')
	rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
				 line_color=None, fill_alpha=0.6)
	p.add_glyph(source, rects)
	p.yaxis.visible = False
	p.grid.visible = False  

	#sequence text view with ability to scroll along x axis
	p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
				x_range=view_range, y_range=ids, tools="xpan,reset",
				min_border=0, toolbar_location='below')#, lod_factor=1)		  
	glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
				text_font="monospace",text_font_size=fontsize)
	rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
				line_color=None, fill_alpha=0.4)
	p1.add_glyph(source, glyph)
	p1.add_glyph(source, rects)

	p1.grid.visible = False
	p1.xaxis.major_label_text_font_style = "bold"
	p1.yaxis.minor_tick_line_width = 0
	p1.yaxis.major_tick_line_width = 0

	p = gridplot([[p],[p1]], toolbar_location='below')
	return p
def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+str(myDict[k])+"\n")
		out.write(k+"\n")
	out.close()
def get_colors(seqs):
	"""make colors for bases in sequence"""
	text = [i for s in list(seqs) for i in s]
	clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white','N':'grey'}
	colors = [clrs[i] for i in text]
	return colors
import uuid
def complex_plot(df,col,expected_sequence,outFigure,width=15,top_alignments=200):
	gRNA_seq = expected_sequence
	aligned_gRNA_seq=""
	fasta_df = pd.DataFrame(df.groupby(col).read_count.sum())
	fasta_df['name'] = fasta_df.index+"_"+fasta_df.read_count.astype(str)
	fasta_df = fasta_df.sort_values("read_count",ascending=False).head(n=top_alignments)
	input = str(uuid.uuid4()).split("-")[-1]
	output = input+".out"
	write_fasta(input,fasta_df.name.to_dict())
	aln = muscle_alignment(input,output) # take a long time
	os.system("rm %s*"%(input))
	# aln  = AlignIO.read("test.out.txt", 'fasta')
	gRNA_list = []
	for i in aln:
		mySplit = i.id.split("_")
		if gRNA_seq == mySplit[0]:
			aligned_gRNA_seq = i.seq
		gRNA_list+=[i.seq]*int(mySplit[-1])
	m = motifs.create([Seq(x.upper()) for x in gRNA_list],alphabet='ACGT-N')
	t=pd.DataFrame.from_dict(m.pwm)
	plt.subplots_adjust(hspace=0)
	fig, (ax1, ax2,ax3) = plt.subplots(3, 1,figsize=(width,6),gridspec_kw={'height_ratios': [1,4, 1]},constrained_layout=True)



	a=sns.heatmap(t.T,annot=True,fmt=".4f",cmap="Blues_r",ax=ax2,cbar=False,xticklabels=False,yticklabels=True)
	ax2.tick_params(axis='y', labelrotation=0 )
	error_rate = []
	for j in range(len(aligned_gRNA_seq)):
		base = aligned_gRNA_seq[j]
		i = "ACGT-N".index(base)
		error = 1-t.at[j,base]
		# print (j,error)
		error_rate.append(error)
		# print (i,j)
		if base!="-":
			ax2.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='red', lw=3))


	ax1.plot(range(1,len(aligned_gRNA_seq)+1),error_rate)	
	ax1.set_xlim(1,len(aligned_gRNA_seq))
	ax1.set_xticks([])
	# print (error_rate)
	logomaker.Logo(t,
				   shade_below=.5,
				   fade_below=.5,
				   color_scheme={'A': 'blue',
					 'C': 'orange',
					 'G': 'green',
					 'T': 'red',
					 '-': 'grey',
					 'N': 'black',
							   },
					ax=ax3)
	# ax3.axis(ymin = 0, ymax = 1, xmin =1, xmax =len(gRNA_seq) )
	ax3.set_xticks([i-1 for i in range(1,len(aligned_gRNA_seq)+1)],[i for i in range(1,len(aligned_gRNA_seq)+1)])
	ax3.set_yticks([])  
	fig.savefig(outFigure,bbox_inches='tight')
	return t

# top 20 alignment plots
def parse_file(f):
	df = pd.read_csv(f,index_col=0)
	df.obs_gRNA = df.obs_gRNA.str.replace("-","")
	df.obs_scaffold = df.obs_scaffold.str.replace("-","")
	df['full_product'] = df.obs_gRNA+df.obs_scaffold
	df['read_count'] = [int(x.split("_")[-1]) for x in df.index]
	df['seq_len'] = df.seq.str.len()	
	return df
files = glob.glob("*alignment.csv")
files
import pickle
from pathlib import Path

# from selenium.webdriver import Firefox
# from selenium.webdriver.firefox.options import Options
# options = Options()
# web_driver = Firefox(options=options,firefox_binary=str(Path('/home/yli11/.conda/envs/captureC/bin/firefox')),executable_path=str(Path("/home/yli11/.conda/envs/captureC/bin/geckodriver")))
# export_png(p, filename= 'hv_plot.png', webdriver=web_driver)

def single_run(f):
	scaffold_seq = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
	df = parse_file(f) 
	print (f)
	gRNA_seq = df.gRNA[0]
	label = f.replace(".alignment.csv","")
	temp_file = "%s.top10.fa"%(label)
	temp_file_out = "%s.top10.out"%(label)
	## gRNA
	tmp = pd.DataFrame(df.groupby("obs_gRNA").read_count.sum().sort_values(ascending=False).head(n=20))
	tmp2 = tmp.read_count.tolist()
	uniquify(tmp2)
	tmp.read_count=tmp2
	
	write_fasta(temp_file,tmp.read_count.to_dict())
	aln = muscle_alignment(temp_file,temp_file_out)
	p=view_alignment(aln, fontsize="9pt", plot_width=500)
	# save(p,)
	# save(p, filename="%s.top20.gRNA.html"%(label))
	export_png(p, filename="%s.top20.gRNA.png"%(label))
	# try:
		# export_png(p, filename="%s.top20.gRNA.png"%(label))
	# except:
		# with open("%s.top20.gRNA.pickle"%(label), 'wb') as handle:
			# pickle.dump(aln, handle, protocol=pickle.HIGHEST_PROTOCOL)

	
	
	## scaffold
	tmp = pd.DataFrame(df.groupby("obs_scaffold").read_count.sum().sort_values(ascending=False).head(n=20))
	tmp2 = tmp.read_count.tolist()
	uniquify(tmp2)
	tmp.read_count=tmp2
	
	write_fasta(temp_file,tmp.read_count.to_dict())
	aln = muscle_alignment(temp_file,temp_file_out)
	p=view_alignment(aln, fontsize="9pt", plot_width=1000)
	# save(p, filename="%s.top20.scaffold.html"%(label))
	export_png(p, filename="%s.top20.scaffold.png"%(label))
	# try:
		# export_png(p, filename="%s.top20.scaffold.png"%(label))
	# except:
		# with open("%s.top20.scaffold.pickle"%(label), 'wb') as handle:
			# pickle.dump(aln, handle, protocol=pickle.HIGHEST_PROTOCOL)
	## full product
	tmp = pd.DataFrame(df.groupby("full_product").read_count.sum().sort_values(ascending=False).head(n=20))
	tmp2 = tmp.read_count.tolist()
	uniquify(tmp2)
	tmp.read_count=tmp2
	
	
	write_fasta(temp_file,tmp.read_count.to_dict())
	aln = muscle_alignment(temp_file,temp_file_out)
	p=view_alignment(aln, fontsize="9pt", plot_width=1200)
	# save(p, filename="%s.top20.full_product.html"%(label))
	export_png(p, filename="%s.top20.full_product.png"%(label))
	# try:
		# export_png(p, filename="%s.top20.full_product.png"%(label))
	# except:
		# with open("%s.top20.full_product.pickle"%(label), 'wb') as handle:
			# pickle.dump(aln, handle, protocol=pickle.HIGHEST_PROTOCOL)


	## read length distribution
	plot_df = pd.DataFrame(df.groupby("seq_len").read_count.sum())
	plt.figure(figsize=(5,20))
	sns.barplot(y=plot_df.index.astype(str).tolist(),x=plot_df.read_count.tolist())
	plt.title("%s read length distribution after trim"%(label))
	plt.savefig("%s.trimmed_read_length_dist.png"%(label),bbox_inches='tight')
	plt.close()
	
	## position error line plot, heatmap overall frequency, sequence logo
	## gRNA
	# df2 = df[(df.gRNA_distance!=-1)&(df.gRNA_distance<=10)&(df.status!="failed")&(df.status!="contained_N")]
	t=complex_plot(df,"obs_gRNA",gRNA_seq,"%s.top200.gRNA.complex_plot.png"%(label),15)
	
	## scaffold
	# df2 = df[(df.scaffold_distance!=-1)&(df.status!="failed")&(df.status!="contained_N")]
	t=complex_plot(df,"obs_scaffold",scaffold_seq,"%s.top200.obs_scaffold.complex_plot.png"%(label),50)
	
	## full product
	# df2 = df[(df.scaffold_distance!=-1)&(df.gRNA_distance!=-1)&(df.status!="failed")&(df.status!="contained_N")]
	t=complex_plot(df,"full_product",gRNA_seq+scaffold_seq,"%s.top200.full_product.complex_plot.png"%(label),65)

def boken_plots(f):
	scaffold_seq = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
	df = parse_file(f) 
	print (f)
	gRNA_seq = df.gRNA[0]
	label = f.replace(".alignment.csv","")
	temp_file = "%s.top10.fa"%(label)
	temp_file_out = "%s.top10.out"%(label)
	## gRNA
	tmp = pd.DataFrame(df.groupby("obs_gRNA").read_count.sum().sort_values(ascending=False).head(n=20))
	tmp2 = tmp.read_count.tolist()
	uniquify(tmp2)
	tmp.read_count=tmp2
	
	write_fasta(temp_file,tmp.read_count.to_dict())
	aln = muscle_alignment(temp_file,temp_file_out)
	p=view_alignment(aln, fontsize="9pt", plot_width=500)
	# save(p,)
	# save(p, filename="%s.top20.gRNA.html"%(label))
	export_png(p, filename="%s.top20.gRNA.png"%(label))
	# try:
		# export_png(p, filename="%s.top20.gRNA.png"%(label))
	# except:
		# with open("%s.top20.gRNA.pickle"%(label), 'wb') as handle:
			# pickle.dump(aln, handle, protocol=pickle.HIGHEST_PROTOCOL)

	
	
	## scaffold
	tmp = pd.DataFrame(df.groupby("obs_scaffold").read_count.sum().sort_values(ascending=False).head(n=20))
	tmp2 = tmp.read_count.tolist()
	uniquify(tmp2)
	tmp.read_count=tmp2
	
	write_fasta(temp_file,tmp.read_count.to_dict())
	aln = muscle_alignment(temp_file,temp_file_out)
	p=view_alignment(aln, fontsize="9pt", plot_width=1000)
	# save(p, filename="%s.top20.scaffold.html"%(label))
	export_png(p, filename="%s.top20.scaffold.png"%(label))
	# try:
		# export_png(p, filename="%s.top20.scaffold.png"%(label))
	# except:
		# with open("%s.top20.scaffold.pickle"%(label), 'wb') as handle:
			# pickle.dump(aln, handle, protocol=pickle.HIGHEST_PROTOCOL)
	## full product
	tmp = pd.DataFrame(df.groupby("full_product").read_count.sum().sort_values(ascending=False).head(n=20))
	tmp2 = tmp.read_count.tolist()
	uniquify(tmp2)
	tmp.read_count=tmp2
	
	
	write_fasta(temp_file,tmp.read_count.to_dict())
	aln = muscle_alignment(temp_file,temp_file_out)
	p=view_alignment(aln, fontsize="9pt", plot_width=1200)
	# save(p, filename="%s.top20.full_product.html"%(label))
	export_png(p, filename="%s.top20.full_product.png"%(label))
	# try:
		# export_png(p, filename="%s.top20.full_product.png"%(label))
	# except:
		# with open("%s.top20.full_product.pickle"%(label), 'wb') as handle:
			# pickle.dump(aln, handle, protocol=pickle.HIGHEST_PROTOCOL)

def plot_complex(f):
	scaffold_seq = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
	df = parse_file(f) 
	print (f)
	gRNA_seq = df.gRNA[0]
	label = f.replace(".alignment.csv","")
	temp_file = "%s.top10.fa"%(label)
	temp_file_out = "%s.top10.out"%(label)

	## read length distribution
	plot_df = pd.DataFrame(df.groupby("seq_len").read_count.sum())
	plt.figure(figsize=(5,20))
	sns.barplot(y=plot_df.index.astype(str).tolist(),x=plot_df.read_count.tolist())
	plt.title("%s read length distribution after trim"%(label))
	plt.savefig("%s.trimmed_read_length_dist.png"%(label),bbox_inches='tight')
	plt.close()
	
	## position error line plot, heatmap overall frequency, sequence logo
	## gRNA
	df2 = df[(df.gRNA_distance!=-1)&(df.gRNA_distance<=10)&(df.status!="failed")&(df.status!="contained_N")]
	t=complex_plot(df2,"obs_gRNA",gRNA_seq,"%s.top200.gRNA.complex_plot.png"%(label),15)
	
	## scaffold
	df2 = df[(df.scaffold_distance!=-1)&(df.status!="failed")&(df.status!="contained_N")]
	t=complex_plot(df2,"obs_scaffold",scaffold_seq,"%s.top200.obs_scaffold.complex_plot.png"%(label),50)
	
	## full product
	df2 = df[(df.scaffold_distance!=-1)&(df.gRNA_distance!=-1)&(df.status!="failed")&(df.status!="contained_N")]
	t=complex_plot(df2,"full_product",gRNA_seq+scaffold_seq,"%s.top200.full_product.complex_plot.png"%(label),65)


# Parallel(n_jobs=4,verbose=10,backend="multiprocessing")(delayed(plot_complex)(i) for i in files)

# for f in files:
	# boken_plots(f)
for f in files:
	single_run(f)
# gRNA mapping summary


def parse_file(f):
	df = pd.read_csv(f,index_col=0)
	df['full_product'] = df.obs_gRNA+df.obs_scaffold
	df['read_count'] = [int(x.split("_")[-1]) for x in df.index]	
	total_read = df.read_count.sum()
	# df = df[df.status!="failed"]
	# df = df[df.status!="containe_N"]
	perfect_scaffold = df[df.scaffold_distance==0].read_count.sum()
	aligned_scaffold = df[(df.scaffold_distance!=-1)&(df.scaffold_distance<=30)].read_count.sum()
	perfect_gRNA = df[df.gRNA_distance==0].read_count.sum()
	aligned_gRNA = df[(df.gRNA_distance!=-1)&(df.gRNA_distance<=10)].read_count.sum()
	perfect_product = df[(df.scaffold_distance==0)&(df.gRNA_distance==0)].read_count.sum()
	aligned_product = df[(df.scaffold_distance!=-1)&(df.scaffold_distance<=30)&(df.gRNA_distance!=-1)&(df.gRNA_distance<=10)].read_count.sum()
	return [total_read,perfect_scaffold,aligned_scaffold,perfect_gRNA,aligned_gRNA,perfect_product,aligned_product]
out = [parse_file(f) for f in files]
stat = pd.DataFrame(out)
stat.columns = ["total_read","perfect_scaffold","aligned_scaffold","perfect_gRNA","aligned_gRNA","perfect_product","aligned_product"]
stat.index = [x.replace(".alignment.csv","") for x in files]
for c in stat.columns:
	if c!="total_read":
		stat[f'{c}_percent'] = stat[c]*100/stat.total_read
stat.to_csv("Mapping_summary.csv")
## scatter plot 

fig, ax = plt.subplots(figsize=(10,10))
sns.scatterplot(data=stat,x="perfect_product_percent",y="perfect_gRNA_percent",ax=ax)
texts = [ax.text(stat.perfect_product_percent.tolist()[i], stat.perfect_gRNA_percent.tolist()[i], stat.index.tolist()[i]) for i in range(stat.shape[0])]
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))
plt.savefig("perfect_gRNA_vs_product_scatter.png",bbox_inches='tight')



# check contaminants
files = glob.glob("*.alignment.csv")
def revcomp(seq):
	try: ## python2
		tab = string.maketrans(b"ACTG", b"TGAC")
	except:  ## python3
		tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]
def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+k+"\n")
		out.write(myDict[k]+"\n")
	out.close()
def check_contaminants(f):
	print (f)
	df = pd.read_csv(f,index_col=0)
	#df.columns = ["gRNA_start","gRNA","obs_gRNA","gRNA_distance","obs_scaffold","scaffold_distance","status","seq"]
	exp_gRNA = df.gRNA.tolist()[0]
	df['g_len'] = df.obs_gRNA.str.len()
	df['read_count'] = [int(x.split("_")[-1]) for x in df.index]
	df['seq_len'] = df.seq.str.len()
	too_short = df[df.seq_len<=25].shape[0]
	df = df[df.seq_len>25]
	tmp = df[df.g_len==len(exp_gRNA)]
	tmp = tmp[(tmp.gRNA_distance<=10)&(tmp.gRNA_distance!=-1)]
	df2 = df[~df.index.isin(tmp.index)]
	df2 = df2.fillna("")
	write_fasta("%s.unaligned.fa"%(f),df2.seq.to_dict())
	df3 = df2.groupby("obs_gRNA")['read_count'].sum().to_dict()
	gRNA_list = ["CTTACCCCACTTAACTATCT","ATTATACCTGCCATGCCGTA","CCCTACCTGTCACCAGGACC","CACTCACCTTAGCCTGAGCA","CACCTACCTAAGAACCATCC","GAACACAAAGCATAGACTGC"]
	for g in ["CTTACCCCACTTAACTATCT","ATTATACCTGCCATGCCGTA","CCCTACCTGTCACCAGGACC","CACTCACCTTAGCCTGAGCA","CACCTACCTAAGAACCATCC","GAACACAAAGCATAGACTGC"]:
		gRNA_list.append(revcomp(g))
	from Levenshtein import distance
	from skbio.alignment import global_pairwise_align_nucleotide as pa
	from skbio import DNA
	scaffold="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
	s2=DNA(scaffold)
	same_gRNA_library_contaminants = 0
	no_gRNA = 0
	unknown = []
	for i in df3:
		o=i
		read_count = df3[i]
		flag = True
		if o=="":
			no_gRNA+=read_count
			continue
		for g in gRNA_list:
			d=distance(g,o)
			if d<=10:
				same_gRNA_library_contaminants+=read_count
				flag = False
				break
		if flag:
			s1 = DNA(o)
			s1_rv = DNA(revcomp(o))
			alignment, score1, start_end_positions = pa(s1,s2)
			alignment, score2, start_end_positions = pa(s1_rv,s2)
			if score1>6 or score2>6:
				no_gRNA+=read_count
			else:
				unknown+=[o]*read_count
	return [f,df2.read_count.sum(),too_short,same_gRNA_library_contaminants,no_gRNA,len(unknown),list(set(unknown))]



out = Parallel(n_jobs=-1,verbose=10,backend="loky")(delayed(check_contaminants)(i) for i in files)
df = pd.DataFrame(out)
casOffinder = [j for sub in df[6].tolist() for j in sub]
casOffinder = pd.DataFrame(list(set(casOffinder)))
print (casOffinder.shape)
casOffinder.to_csv("unaligned.casOffinder_to_check.list",index=False,header=False)
df = df.set_index(0)
df.columns = ['total_unaligned_gRNA_reads','too_short','same_gRNA_library_contaminants_counts','no_gRNA_counts','for_casOffinder_counts','for_casOffinder']
df.to_csv("unaligned.stat.csv")
