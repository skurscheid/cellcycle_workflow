__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-015"

from snakemake.exceptions import MissingInputException
import os

REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"]
RUN_ID = "N08851_SK_LR1807201_SEQ"
PROJECT_ID = "LR1807201"
ASSAY_TYPE = "ChIP-Seq"

#singularity: "docker://skurscheid/snakemake_baseimage:0.1"
#singularity: "docker://continuumio/miniconda3:4.4.10"

rule:
    version:
        "1.0"

localrules:
        all

home = os.environ['HOME']

include_prefix = "rules/"

include:
      include_prefix + "run_macs2.smk"

rule callpeaks:
    input:
        expand("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-{rep}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["G1"],
                treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
		rep = ["1", "2"]),
        expand("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-{rep}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["M"],
                treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                rep = ["1", "2"])
