#random_seed = config.get("seed", 2020)

DATA = config.get("data",None)

TOPGO_NODES = config.get("topgo_nodes",100)

MIN_QVAL = int(config.get("min_qval",0.005) * 1000)
MAX_QVAL = int(config.get("max_qval",0.1) * 1000)
STEP_QVAL = int(config.get("step_qval", 0.005) * 1000)

MIN_COMPS = config.get("min_comps",50)
MAX_COMPS = config.get("max_comps",500) + 1
STEP_COMPS = config.get("step_comps", 10)

REPS = range(1,config.get("reps",1) + 1,1)

COMPONENTS = range(MIN_COMPS,MAX_COMPS,STEP_COMPS)
QVALS = [x / 1000.0 for x in range(MIN_QVAL,MAX_QVAL, STEP_QVAL)]

if config.get("is_test",False):
    REPS = [1]
    COMPONENTS = [20]
    QVALS = [0.01]

RELATIONS = ["spearman-abs","pearson-abs","bicor-tom-abs",
             "spearman-signed","pearson-signed","bicor-tom-signed"]

RELATIONS = ["bicor-tom-abs","pearson-abs",
             "bicor-tom-signed","pearson-signed"]


ONTS = ["BP","MF","CC"]

ICA_VERSIONS = [1,2]

rule target:
    input:
        expand("ica_fdr{v}_{c}comps_rep{rep}_qvalues.csv.gz",c=COMPONENTS,rep=REPS,v=ICA_VERSIONS),
        expand("ica_fdr{v}_{c}comps_rep{rep}_{q}qval.json", c=COMPONENTS, q=QVALS, rep=REPS,v=ICA_VERSIONS),
        expand("enr_fdr{v}_{c}comps_rep{r}_{f}qval_{o}.csv",c=COMPONENTS,r=REPS,f=QVALS,o=ONTS,v=ICA_VERSIONS),
        "enr.pdf",
        "igp.pdf"

# ------------------------------------------------------------------------------
# pre-processing
# ------------------------------------------------------------------------------

rule standardize:
    input:
        DATA
    output:
        "standardized.feather"
    conda:
        "envs/all.yaml"
    script:
        "scripts/standardize.py"

# ------------------------------------------------------------------------------
# calculate correlations
# ------------------------------------------------------------------------------


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

# ------------------------------------------------------------------------------
# ICA itself
# ------------------------------------------------------------------------------

rule ica:
    input:
        rules.standardize.output,
    output:
        "ica_{components}comps_rep{rep}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/ica.py"

# ------------------------------------------------------------------------------
# fdr calcs and cutoffs
# ------------------------------------------------------------------------------

rule fdr_calc:
    input:
        rules.ica.output,
    output:
        "ica_fdr{ICAver}_{components}comps_rep{rep}_qvalues.csv.gz"
    params:
        ICAver = lambda wc: wc.ICAver
    conda:
        "envs/all.yaml"
    script:
        "scripts/qval_calc.R"

rule fdr_cut:
    input:
        rules.fdr_calc.output,
    output:
        "ica_fdr{ICAver}_{components}comps_rep{rep}_{fdr}qval.json"
    params:
        q = lambda wc: float(wc.fdr)
    conda:
        "envs/all.yaml"
    script:
        "scripts/qval_cutoff_1d.R"

# ------------------------------------------------------------------------------
# in-group proportion metrics
# ------------------------------------------------------------------------------

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
        "igp_fdr{ICAver}_{components}comps_rep{rep}_{fdr}qval_{relation}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/igp.R"

rule plot_igp_maximization:
    input:
        lambda wc: expand("igp_fdr{v}_{c}comps_rep{rep}_{q}qval_{r}.csv", c=COMPONENTS, q=QVALS, r=RELATIONS, rep=REPS, v=ICA_VERSIONS),
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
        "enr_fdr{ICAver}_{components}comps_rep{rep}_{fdr}qval_{ont}.csv"
    params:
        qval = lambda wc: wc.fdr,
        nodes = TOPGO_NODES,
    #conda:
    #    "envs/topgo.yaml"
    script:
        "scripts/go_enr.R"

rule plot_enr_maximization:
    input:
        expand("enr_fdr{v}_{c}comps_rep{r}_{f}qval_{o}.csv",c=COMPONENTS,r=REPS,f=QVALS,o=ONTS, v=ICA_VERSIONS)
    output:
        "enr.pdf"
    conda:
        "envs/all.yaml"
    script:
        "scripts/plot_enr.R"

rule mgc:
    input:
        communities = rules.fdr_cut.output,
        corrmat = "{relation}.feather",
    output:
        "mgc_fdr{ICAver}_{components}comps_rep{rep}_{fdr}qval_{relation}.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/mgc.R"

#rule plot_mgc_maximization:
#    input:
#        expand("mgc_fdr{v}_{c}comps_rep{rep}_{q}qval_{r}.csv", c=COMPONENTS, q=QVALS, r=RELATIONS, rep=REPS, v=ICA_VERSIONS),
#    output:
#        "mgc.pdf"
#    conda:
#        "envs/all.yaml"
#    script:
#        "scripts/plot_mgc.R"
