

import pandas as pd
import sys
import os
import glob
aligned_csv = sys.argv[1]
label=sys.argv[2]
try:
	df = pd.read_csv(aligned_csv,index_col=0)
except:
	aligned_csv = glob.glob(aligned_csv)[0]
	df = pd.read_csv(aligned_csv,index_col=0)
df = df[df.aligned_gRNA.isnull()]
df["obs_gRNA_length"] = df.obs_gRNA.str.len()
df = df[df.obs_gRNA_length.between(17,23)]
df['read_count'] = [int(x.split("_")[-1]) for x in df.index]
df2 = df.groupby('obs_gRNA')['read_count'].sum().sort_values().reset_index()
df2.columns = ['gRNA_noise','counts']
df2 = df2[df2.counts>2]
# df2 = df.obs_gRNA.value_counts().reset_index()
# df2 = df2[df2.obs_gRNA>2]
# df2.columns = ['gRNA_noise','counts']
# from joblib import Parallel, delayed
print (df2)
exit()
for i,r in df2.iterrows():
	print (r)
	command = f'calitas SearchReference -i {r.gRNA_noise}ngg -I {r.gRNA_noise}_{r.counts} -r /home/yli11/Data/Human/hg38/fasta/hg38.main.fa -o {label}_{r.gRNA_noise}.csv --max-guide-diffs 0 --max-pam-mismatches 0 --max-gaps-between-guide-and-pam 0'
	os.system(command)


# Parallel(n_jobs=-1,verbose=10)(delayed(os.system)(f'calitas SearchReference -i {r.gRNA_noise}ngg -I {r.gRNA_noise}_{r.counts} -r /home/yli11/Data/Human/hg38/fasta/hg38.main.fa -o {label}_{r.gRNA_noise}.csv --max-guide-diffs 0 --max-pam-mismatches 0 --max-gaps-between-guide-and-pam 0') for r in df2.to_dict(orient="records"))


