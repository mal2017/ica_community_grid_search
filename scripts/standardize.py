import pandas as pd
import numpy as np

def standardize(X, maxval):
    Z = (X - X.mean())/(X.std())
    if maxval:
      return Z.clip(-1*maxval,maxval)
    else:
      return Z

expression = pd.read_feather(snakemake.input[0])
expression.set_index('index', inplace=True)
expression = standardize(expression, snakemake.params["maxval"]).T
expression.reset_index().to_feather(snakemake.output[0])
