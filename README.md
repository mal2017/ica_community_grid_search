
## Data

Feed depth normalized, logarithmized matrix in csv format. Genes x samples.

```
--config data=some/matrix.csv.gz
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
		--config data=some/matrix.csv.gz
```

## Testing

```
--config is_test=True
```
