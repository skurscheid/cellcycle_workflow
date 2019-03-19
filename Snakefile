__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-015"

from snakemake.exceptions import MissingInputException
import os
from snakemake.remote.AzureStorage import RemoteProvider as AzureRemoteProvider

# setup Azure Storage for remote access
account_key=os.environ['AZURE_KEY']
account_name=os.environ['AZURE_ACCOUNT']
AS = AzureRemoteProvider(account_name=account_name, account_key=account_key)


REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"]
RUN_ID = "N08851_SK_LR1807201_SEQ"
PROJECT_ID = "LR1807201"
ASSAY_TYPE = "ChIP-Seq"

#singularity: "docker://skurscheid/snakemake_baseimage:0.1"
#singularity: "docker://continuumio/miniconda3:4.4.10"

wildcard_constraints:
    chip_library = ".*?(?=\_)"

rule:
    version:
        "1.0"

localrules:
       download_treatment, download_control

home = os.environ['HOME']

include_prefix = "rules/"

include:
      include_prefix + "deepTools_data_prep.smk"
include:
      include_prefix + "azure_rules.smk"
include:
      include_prefix + "run_macs2.smk"
include:
      include_prefix + "run_idr.smk"

rule callpeaks:
    input:
        expand("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["G1"],
                treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                suffix = ["1", "2"]),
        expand("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["M"],
                treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                suffix = ["1", "2"])

rule idr:
    input:
        expand("{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{treatment}/{treatment}_{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["G1"],
                treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                suffix = ["idr.bed", "other.bed"]),
        expand("{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{treatment}/{treatment}_{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["M"],
                treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                suffix = ["idr.bed", "other.bed"])

rule plot_peaks_per_sample_example:
    input:
        AS.remote( expand("experiment/{assayType}/{project}/{runID}/deepTools/plotHeatmap/{reference_version}/{cycle}/{library}-{rep}.pdf",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["M"],
                library = "ACTR6M",
                rep = ["1", "2"],
                suffix = ["pdf"]) )

rule peaks_per_sample:
        input:
            expand("{assayType}/{project}/{runID}/deepTools/plotHeatmap/{reference_version}/{cycle}/{treatment}-{rep}.pdf",
                   assayType = "ChIP-Seq",
                   reference_version = REF_VERSION,
                   project = "LR1807201",
                   runID = "N08851_SK_LR1807201_SEQ",
                   cycle = ["G1"],
                   treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                   rep = ["1", "2"],
                   suffix = ["pdf"]),
            expand("{assayType}/{project}/{runID}/deepTools/plotHeatmap/{reference_version}/{cycle}/{treatment}-{rep}.pdf",
                   assayType = "ChIP-Seq",
                   reference_version = REF_VERSION,
                   project = "LR1807201",
                   runID = "N08851_SK_LR1807201_SEQ",
                   cycle = ["M"],
                   treatment = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                   rep = ["1", "2"],
                   suffix = ["pdf"])
