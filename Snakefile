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

#RELATIONS = ["bicor-tom-abs"]

FUNCTIONAL = ["top5Enr"]

ONTS = ["BP","MF","CC"]

rule target:
    input:
        #expand("neighbors_{r}.json",r=RELATIONS),
        #expand("ica_{c}comps_rep{rep}_{q}qval.json", c=COMPONENTS, q=QVALS, rep=REPS),
        #expand("ica_{c}comps_rep{rep}_qvalues.csv.gz",c=COMPONENTS,rep=REPS),
        "igp.pdf",
        expand("enr_{c}comps_rep{r}_{f}qval_{o}.csv",c=COMPONENTS,r=REPS,f=QVALS,o=ONTS),
        "enr.pdf"

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

rule fdr_calc:
    input:
        rules.ica.output,
    output:
        "ica_{components}comps_rep{rep}_qvalues.csv.gz"
    conda:
        "envs/all.yaml"
    script:
        "scripts/qval_calc.R"

rule fdr_cut:
    input:
        rules.fdr_calc.output,
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
        communities = rules.fdr_cut.output,
        nn = rules.neighbors.output,
    output:
        "igp_{components}comps_rep{rep}_{fdr}qval_{relation}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/igp.R"

rule plot_igp_maximization:
    input:
        lambda wc: expand("igp_{c}comps_rep{rep}_{q}qval_{r}.csv", c=COMPONENTS, q=QVALS, r=RELATIONS, rep=REPS),
    output:
        "igp.pdf"
    conda:
        "envs/all.yaml"
    script:
        "scripts/plot_igp.R"

rule run_topgo:
    input:
        rules.fdr_calc.output
    output:
        "enr_{components}comps_rep{rep}_{fdr}qval_{ont}.csv"
    params:
        qval = lambda wc: wc.fdr
    #conda:
    #    "envs/topgo.yaml"
    script:
        "scripts/go_enr.R"

rule plot_enr_maximization:
    input:
        expand("enr_{c}comps_rep{r}_{f}qval_{o}.csv",c=COMPONENTS,r=REPS,f=QVALS,o=ONTS)
    output:
        "enr.pdf"
    conda:
        "envs/all.yaml"
    script:
        "scripts/plot_enr.R"
