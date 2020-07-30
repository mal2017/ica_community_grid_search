import pandas as pd
import numpy as np
import random
from sklearn.decomposition import FastICA

seed = snakemake.params['random_seed']
random.seed(seed)
np.random.seed(seed)

X = pd.read_csv(snakemake.input[0])
X.set_index('index',inplace=True)
ica = FastICA(n_components=int(snakemake.params['comps']), random_state=seed, max_iter=1000)
source = ica.fit_transform(X)

mixing = pd.DataFrame(ica.mixing_, X.columns)
res = pd.DataFrame(source, X.index)

mixing.to_csv(snakemake.output["mixing"])
res.to_csv(snakemake.output["source"])
