__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-015"

from snakemake.exceptions import MissingInputException
import os

configfile: "config.json"

REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"]
RUN_ID = "N08851_SK_LR1807201_SEQ"
PROJECT_ID = "LR1807201"
dataDir = "/data/cellcycledata/"

#singularity: "docker://skurscheid/snakemake_baseimage:0.1"
singularity: "docker://continuumio/miniconda3:4.4.10"

rule:
    version:
        "1.0"

localrules:
    all

home = os.environ['HOME']

include_prefix = "rules/"

include:
    include_prefix + "run_fastp.smk"
include:
    include_prefix + "run_alignment.smk"
include:
    include_prefix + "deepTools_QC.smk"
include:
    include_prefix + "deepTools_plotting.smk"

rule execute_collectInsertSize:
    input:
        expand("{assayType}/picardTools/CollectInsertSizeMetrics/{reference_version}/{runID}/{library}.{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                library = [x for x in config["samples"]["ChIP-Seq"]["LR1807201"]["N08851_SK_LR1807201_SEQ"].keys()],
                suffix = ["histogram.pdf", "insert_size_metrics.txt"])

rule execute_fastp:
    input:
        expand(dataDir + "{assayType}/{project}/{runID}/fastp/trimmed/{library}.{suffix}",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               reference_version = REF_VERSION,
               runID = RUN_ID,
               library = [x for x in config["samples"]["ChIP-Seq"]["LR1807201"]["N08851_SK_LR1807201_SEQ"].keys()],
	           suffix = ["end1.fastq.gz", "end2.fastq.gz"]),
        expand(dataDir + "{assayType}/{project}/{runID}/fastp/report/{library}.{suffix}",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               reference_version = REF_VERSION,
               runID = RUN_ID,
               library = [x for x in config["samples"]["ChIP-Seq"]["LR1807201"]["N08851_SK_LR1807201_SEQ"].keys()],
	           suffix = ["fastp.html", "fastp.json"])

rule execute_deepTools_QC:
    input:
        expand("{assayType}/{project}/{runID}/deepTools/plotFingerprint/{reference_version}/{condition}/fingerprints.png",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               runID = RUN_ID,
               reference_version = REF_VERSION,
               condition = ["G1", "M"]),
        expand("{assayType}/{project}/{runID}/deepTools/bamPEFragmentSize/{reference_version}/{condition}/histogram.png",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               runID = RUN_ID,
               reference_version = REF_VERSION,
               condition = ["G1", "M"]),
        expand("{assayType}/{project}/{runID}/deepTools/plotPCA/{reference_version}/{condition}/PCA_readCounts.png",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               runID = RUN_ID,
               reference_version = REF_VERSION,
               condition = ["G1", "M"]),
        expand("{assayType}/{project}/{runID}/deepTools/plotCorrelation/{reference_version}/{condition}/heatmap_SpearmanCorr_readCounts.{suffix}",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               runID = RUN_ID,
               reference_version = REF_VERSION,
               condition = ["G1", "M"],
               suffix = ["png", "tab"])


rule execute_deepTools_plotting:
    input:
        expand("{assayType}/{project}/{runID}/deepTools/plotProfile/{subcommand}/{reference_version}/{region}_{suffix}.pdf",
               assayType = "ChIP-Seq",
               project = PROJECT_ID,
               reference_version = REF_VERSION,
               runID = RUN_ID,
               subcommand = "scale-region",
               region = ["allGenes"],
	       suffix = "RPKM")


rule all:
    input:
        expand("{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{library}.bam.bai",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                library = [x for x in config["samples"]["ChIP-Seq"]["LR1807201"]["N08851_SK_LR1807201_SEQ"].keys()])
