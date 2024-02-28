
import scanpy as sc
from scipy import io
import sys
import os
import re 
from scipy.sparse import csc_matrix
import gzip 

os.chdir('covid_szabo/sims-farber/')
os.getcwd()
files = os.listdir()
r = re.compile(".*h5ad")
file_h5ad = list(filter(r.match, files))

for sample in file_h5ad:
  out_dir =  sample.split('.', 1)[0]
  print(out_dir)
  adata = sc.read_h5ad(sample)
  os.makedirs(out_dir,exist_ok=True)
  adata = adata.raw.to_adata() #only if adata has RAW saved and thats what you want!!

  with open(out_dir + '/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
      f.write(item + '\n')

  with open(out_dir + '/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
      f.write(item + '\n')
    
  io.mmwrite(out_dir +'/matrix', csc_matrix(adata.X.T))

  os.system("gzip -r "+out_dir + '/*')

  adata.obs.to_csv(out_dir + '/metadata.csv')
