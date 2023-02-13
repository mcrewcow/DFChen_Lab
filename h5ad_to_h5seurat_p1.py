import scanpy as sc
import numpy as np
import pandas as pd

print(sc.__version__)
print(np.__version__)
print(pd.__version__)

print("Reading data...")
adata = sc.read_h5ad("file.h5ad")

print("Writing data...")
t=adata.raw.X.toarray()
pd.DataFrame(data=t, index=adata.obs_names, columns=adata.raw.var_names).to_csv('raw-counts.csv')
# comma separated csv with header, no rownames
pd.DataFrame(adata.obs).to_csv("metadata.csv")
