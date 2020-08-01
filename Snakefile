#container: "docker://bioconductor/bioconductor_docker:RELEASE_3_11"
import random

RANDOM_SEED = config.get("RANDOM_SEED",2020)
random.seed(RANDOM_SEED)

DATA = config.get("data",None)

TOPGO_NODES = config.get("topgo_nodes",100)

MIN_QVAL = int(config.get("min_qval",0.005) * 1000)
MAX_QVAL = int(config.get("max_qval",0.1) * 1000)
STEP_QVAL = int(config.get("step_qval", 0.01) * 1000)

MIN_COMPS = config.get("min_comps",20)
MAX_COMPS = config.get("max_comps",300) + 1
STEP_COMPS = config.get("step_comps", 10)

REPS = range(1,config.get("reps",10) + 1,1)
OVERALL_REPS = range(1,config.get("overall_reps",3) + 1,1)

COMPONENTS = range(MIN_COMPS,MAX_COMPS,STEP_COMPS)
QVALS = [x / 1000.0 for x in range(MIN_QVAL,MAX_QVAL, STEP_QVAL)]

if config.get("is_test",False):
    REPS = [1,2]
    OVERALL_REPS = [1,2]
    COMPONENTS = [12,14,16,18]
    QVALS = [0.005,0.01,0.05]

INDIV_RANDOM_SEEDS_ICA = [random.sample(range(0,100000,1),k=max(REPS)) for x in OVERALL_REPS]
OUTLIER_FILT_KNN_ICA = config.get("OUTLIER_FILT_KNN_ICA",1)
OUTLIER_MAX_DIST_ICA = config.get("OUTLIER_MAX_DIST_ICA",1300)

ONT = "CC"

ICA_VERSION = 1

print(INDIV_RANDOM_SEEDS_ICA)

rule target:
    input:
        #expand("ica/{c}/{ovr}/{r}/source.csv.gz",c=COMPONENTS, ovr = OVERALL_REPS, r=REPS),
        #expand("ica/{c}/{ovr}/{r}/mixing.csv.gz",c=COMPONENTS, ovr = OVERALL_REPS, r=REPS),
        #expand("ica/{c}/{ovr}/consensus-ica.csv.gz", c=COMPONENTS, ovr=OVERALL_REPS),
        #expand("ica/{c}/{ovr}/consensus-usage.csv.gz", c=COMPONENTS, ovr=OVERALL_REPS),
        expand("ica/{c}/{ovr}/consensus-ica.qvalues.csv.gz",c=COMPONENTS, ovr=OVERALL_REPS),
        expand("enr/{c}/{ovr}/{fdr}/enr.csv",c=COMPONENTS, ovr=OVERALL_REPS, fdr = QVALS),
        "enr/enr.csv",
        "enr/enr-metrics.csv",
        "ica/silhouette.csv",
        #"metrics.csv"

# ------------------------------------------------------------------------------
# ICA itself
# ------------------------------------------------------------------------------

rule ica_reps:
    input:
        DATA
    output:
        source = "ica/{k}/{ovr}/{rep}/source.csv.gz",
        mixing = "ica/{k}/{ovr}/{rep}/mixing.csv.gz",
    params:
        random_seed = lambda wc: INDIV_RANDOM_SEEDS_ICA[int(wc.ovr)-1][int(wc.rep)-1],
        comps = lambda wc: int(wc.k)
    conda:
        "envs/all.yaml"
    script:
        "scripts/ica.py"

rule ica_consensus:
    input:
        source= lambda wc: expand("ica/{k}/{ovr}/{rep}/source.csv.gz",k=wc.k,ovr = wc.ovr, rep=range(1,max(REPS)+1,1)),
        mixing= lambda wc: expand("ica/{k}/{ovr}/{rep}/mixing.csv.gz",k=wc.k, ovr = wc.ovr, rep=range(1,max(REPS)+1,1)),
    output:
        usage="ica/{k}/{ovr}/consensus-usage.csv.gz",
        ica="ica/{k}/{ovr}/consensus-ica.csv.gz",
        dists = "ica/{k}/{ovr}/consensus-dists.csv.gz",
        silhouette = "ica/{k}/{ovr}/consensus-silhouette.csv.gz",
    params:
        knn = OUTLIER_FILT_KNN_ICA,
        max_dist = OUTLIER_MAX_DIST_ICA,
        k = lambda wc: wc.k
    conda:
        "envs/all.yaml"
    script:
        "scripts/ica-consensus.R"

# # ------------------------------------------------------------------------------
# # fdr calcs and cutoffs
# # ------------------------------------------------------------------------------
#
rule fdr_calc:
    input:
        rules.ica_consensus.output.ica,
    output:
        "ica/{k}/{ovr}/consensus-ica.qvalues.csv.gz"
    params:
        ICAver = ICA_VERSION
    conda:
        "envs/all.yaml"
    script:
        "scripts/qval_calc.R"

# # ------------------------------------------------------------------------------
# # enrichment
# # ------------------------------------------------------------------------------

rule run_topgo:
    input:
        rules.fdr_calc.output
    output:
        "enr/{k}/{ovr}/{fdr}/enr.csv"
    params:
        qval = lambda wc: wc.fdr,
        nodes = TOPGO_NODES,
        ont = ONT
    #conda:
        #"envs/topgo.yaml"
    script:
        "scripts/go_enr.R"

rule combine_reps_enr:
    input:
        expand("enr/{k}/{ovr}/{fdr}/enr.csv",k=COMPONENTS,ovr=OVERALL_REPS,fdr=QVALS)
    output:
        csv="enr/enr.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/plot_enr.R"

rule combine_enr_metrics:
    input:
      rules.combine_reps_enr.output
    output:
      "enr/enr-metrics.csv"
    #conda:
        #"envs/topgo.yaml"
    script:
      "scripts/gather_supp_enr_metrics.R"

rule combine_silhouette:
    input:
        expand("ica/{k}/{ovr}/consensus-silhouette.csv.gz",k=COMPONENTS,ovr=OVERALL_REPS)
    output:
        sil = "ica/silhouette.csv"
    conda:
        "envs/all.yaml"
    script:
        "scripts/gather_silhouette.R"

rule reconstruct:
