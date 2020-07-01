
## Data

Feed depth normalized, logarithmized matrix in feather format. Samples x genes.

```
--config data=some/matrix.feather
```

## Get started

```
snakemake --use-conda \
	  --cores 999 \  
	  -u amarel.json \
	  --cluster "sbatch --export=ALL --partition {cluster.partition} \
			--nodes {cluster.n} --time {cluster.time} --ntasks {cluster.tasks} \
		 	--mem {cluster.mem} \
	  	--cpus-per-task={cluster.cpus}" \
		-d ~/scratch/mal456/[something] \
		--config data=some/matrix.feather
```

## Testing

```
--config is_test=True
```
