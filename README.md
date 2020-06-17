snakemake --use-conda \
	  --cores 8 \  
	  -u amarel.json \
	  --cluster "sbatch --export=ALL --partition {cluster.partition} --nodes {cluster.n} --time {cluster.time} --ntasks {cluster.tasks} --mem {cluster.mem} \
	   --cpus-per-task={cluster.cpus}" \
		 -d ~/scratch/mal456/[something] \
		 -kp
