random_seed = config.get("seed", 2020)

DATA = config.get("data",None)

MIN_QVAL = int(config.get("min_qval",0.002) * 1000)
MAX_QVAL = int(config.get("max_qval",0.02) * 1000)
STEP_QVAL = int(config.get("step_qval", 0.002) * 1000)

MIN_COMPS = config.get("min_comps",2)
MAX_COMPS = config.get("max_comps",10)
STEP_COMPS = config.get("step_comps", 2)

COMPONENTS = range(MIN_COMPS,MAX_COMPS,STEP_COMPS)
QVALS = [x / 1000.0 for x in range(MIN_QVAL,MAX_QVAL, STEP_QVAL)]

rule target:
    input:
        expand("neighbors_{r}.json",r=["spearman","pearson"]),
        expand("ica_{c}comps_{q}qval.json", c=COMPONENTS, q=QVALS),
        expand("igp_{r}_grid.csv",r=["spearman","pearson"]),
        "standardized.feather",
        expand("ica_{c}comps.csv",c = COMPONENTS),

rule standardize:
    input:
        DATA
    output:
        "standardized.feather"
    conda:
        "envs/all.yaml"
    script:
        "scripts/standardize.py"

rule pearson:
    input:
        rules.standardize.output,
    output:
        "pearson.feather"
    conda:
        "envs/all.yaml"
    script:
        "scripts/pearson.py"

rule spearman:
    input:
        rules.standardize.output,
    output:
        "spearman.feather"
    conda:
        "envs/all.yaml"
    script:
        "scripts/spearman.py"


rule ica:
    input:
        rules.standardize.output,
    output:
        "ica_{components}comps.csv"
    params:
        random_seed = random_seed
    conda:
        "envs/all.yaml"
    script:
        "scripts/ica.py"

rule fdr:
    input:
        rules.ica.output,
    output:
        "ica_{components}comps_{fdr}qval.json"
    params:
        seed = random_seed,
        q = lambda wc: float(wc.fdr)
    conda:
        "envs/all.yaml"
    script:
        "scripts/qval_cutoff_1d.R"

rule neighbors:
    input:
        "{relation}.feather"
    output:
        "neighbors_{relation}.json"
    conda:
        "envs/all.yaml"
    script:
        "scripts/neighbors_correlation.R"

rule igp:
    input:
        communities = rules.fdr.output,
        nn = rules.neighbors.output,
    output:
        "igp_{components}comps_{fdr}qval_{relation}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/igp.R"

rule collect_igp:
    input:
        lambda wc: expand("igp_{c}comps_{q}qval_{r}.csv", c=COMPONENTS, q=QVALS, r=wc.relation),
    output:
        "igp_{relation}_grid.csv"
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -q -n +2 {input} >> {output}
        """
