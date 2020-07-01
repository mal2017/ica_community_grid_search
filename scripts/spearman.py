import pandas as pd
import numpy as np

X = pd.read_feather(snakemake.input[0])
X.set_index('index', inplace=True)
X = X.T

X = X.rank()
C = X.corr()
Cabs = C.abs()

C.reset_index().to_feather(snakemake.output["signed"])
Cabs.reset_index().to_feather(snakemake.output["abs"])
