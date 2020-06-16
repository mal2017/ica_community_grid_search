import pandas as pd
import numpy as np

def standardize(X):
    return (X - X.mean())/(X.std())

expression = pd.read_feather(snakemake.input[0])
expression.set_index('index', inplace=True)
expression = standardize(expression).T
expression.reset_index().to_feather(snakemake.output[0])
