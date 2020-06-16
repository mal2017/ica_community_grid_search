import pandas as pd
import numpy as np

X = pd.read_feather(snakemake.input[0])
X.set_index('index', inplace=True)
C = X.corr()
C.reset_index().to_feather(snakemake.output[0])
