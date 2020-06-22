#random_seed = config.get("seed", 2020)

DATA = config.get("data",None)

MIN_QVAL = int(config.get("min_qval",0.005) * 1000)
MAX_QVAL = int(config.get("max_qval",0.1) * 1000)
STEP_QVAL = int(config.get("step_qval", 0.005) * 1000)

MIN_COMPS = config.get("min_comps",50)
MAX_COMPS = config.get("max_comps",500) + 1
STEP_COMPS = config.get("step_comps", 10)

REPS = range(1,config.get("reps",1) + 1,1)

COMPONENTS = range(MIN_COMPS,MAX_COMPS,STEP_COMPS)
QVALS = [x / 1000.0 for x in range(MIN_QVAL,MAX_QVAL, STEP_QVAL)]

RELATIONS = ["spearman-abs","pearson-abs","bicor-tom-abs",
             "spearman-signed","pearson-signed","bicor-tom-signed"]


#print(expand("ica_{c}comps_rep{rep}_{q}qval.json", c=COMPONENTS, q=QVALS, rep=REPS))
#print(QVALS)

rule target:
    input:
        expand("neighbors_{r}.json",r=RELATIONS),
        expand("ica_{c}comps_rep{rep}_{q}qval.json", c=COMPONENTS, q=QVALS, rep=REPS),
        expand("igp_{r}_grid.csv",r=RELATIONS),
        "standardized.feather",
        expand("ica_{c}comps_rep{rep}.csv",c = COMPONENTS,rep=REPS),
        "igp.pdf"

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
        signed = "pearson-signed.feather",
        abs = "pearson-abs.feather"
    conda:
        "envs/all.yaml"
    script:
        "scripts/pearson.py"

rule spearman:
    input:
        rules.standardize.output,
    output:
        signed="spearman-signed.feather",
        abs="spearman-abs.feather"
    conda:
        "envs/all.yaml"
    script:
        "scripts/spearman.py"

rule bicor_tom:
    input:
        rules.standardize.output,
    output:
        signed = "bicor-tom-signed.feather",
        abs = "bicor-tom-abs.feather"
    conda:
        "envs/all.yaml"
    threads:
        4
    script:
        "scripts/bicor_tom.R"

rule ica:
    input:
        rules.standardize.output,
    output:
        "ica_{components}comps_rep{rep}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/ica.py"

rule fdr:
    input:
        rules.ica.output,
    output:
        "ica_{components}comps_rep{rep}_{fdr}qval.json"
    params:
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
        "igp_{components}comps_rep{rep}_{fdr}qval_{relation}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/igp.R"

rule collect_igp:
    input:
        lambda wc: expand("igp_{c}comps_rep{rep}_{q}qval_{r}.csv", c=COMPONENTS, q=QVALS, r=wc.relation, rep=REPS),
    output:
        "igp_{relation}_grid.csv"
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -q -n +2 {input} >> {output}
        """

rule plot_igp_maximization:
    input:
        expand("igp_{r}_grid.csv",r=RELATIONS),
    output:
        "igp.pdf"
    conda:
        "envs/all.yaml"
    script:
        "scripts/plot_igp.R"
