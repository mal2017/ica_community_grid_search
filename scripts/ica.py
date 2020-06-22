import pandas as pd
import numpy as np
import random
from sklearn.decomposition import FastICA

#random.seed(snakemake.params.random_seed)
#np.random.seed(snakemake.params.random_seed)

X = pd.read_feather(snakemake.input[0])
X.set_index('index',inplace=True)
X = X.T
ica = FastICA(n_components=int(snakemake.wildcards.components))
source = ica.fit_transform(X)
res = pd.DataFrame(source, X.index)
res.to_csv(snakemake.output[0])
