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
       download_chip_library, download_control

home = os.environ['HOME']

include_prefix = "./rules/"

include:
      include_prefix + "azure_rules.smk"
include:
      include_prefix + "run_macs2.smk"
include:
      include_prefix + "run_idr.smk"

rule callpeaks:
    input:
        expand("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{chip_library}-{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["G1"],
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                suffix = ["1", "2"]),
        expand("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{chip_library}-{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["M"],
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                suffix = ["1", "2"])

rule idr:
    input:
        expand("{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{chip_library}/{chip_library}_{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["G1"],
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                suffix = ["idr.bed", "other.bed"]),
        expand("{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{chip_library}/{chip_library}_{suffix}",
                assayType = "ChIP-Seq",
                reference_version = REF_VERSION,
                project = "LR1807201",
                runID = "N08851_SK_LR1807201_SEQ",
                cycle = ["M"],
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                suffix = ["idr.bed", "other.bed"])

rule peaks_per_sample:
        input:
            expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/{reference_version}/G1/{chip_library}-{rep}.pdf",
                   reference_version = REF_VERSION,
                   chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                   rep = ["1", "2"]),
            expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/{reference_version}/M/{chip_library}-{rep}.pdf",
                   reference_version = REF_VERSION,
                   chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                   rep = ["1", "2"])

rule plot_peaks_per_sample_example:
    input:
        AS.remote("experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/GRCh38_ensembl84/G1/ANP32G1-1.pdf")
