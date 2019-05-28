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

wildcard_constraints:
    chip_library = ".*?(?=\_)"

rule:
    version:
        "1.0"

localrules:
       all

home = os.environ['HOME']

# include external python functions
include:
        "./scripts/common/input_functions.py"
include:
        "./scripts/common/parameter_functions.py"

rule all_callpeaks:
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

rule all_idr:
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

rule all_peaks_per_sample:
    input:
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/{reference_version}/G1/{chip_library}-{rep}.pdf",
                reference_version = REF_VERSION,
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                rep = ["1", "2"]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/{reference_version}/M/{chip_library}-{rep}.pdf",
                reference_version = REF_VERSION,
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                rep = ["1", "2"])

rule all_multimap:
    input:
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/picardTools/MarkDuplicates/{chip_library}-{rep}_mm.bam",
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                rep = ["1", "2"]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/picardTools/MarkDuplicates/{chip_library}-{rep}_mm.bam",
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                rep = ["1", "2"]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/picardTools/MarkDuplicates/{input_library}-{rep}_mm.bam",
                input_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["Input"].keys()],
                rep = ["1", "2"]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/picardTools/MarkDuplicates/{input_library}-{rep}_mm.bam",
                input_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["Input"].keys()],
                rep = ["1", "2"])


rule test_multimap:
    input:
        "ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/picardTools/MarkDuplicates/H2AZM-1_mm.bam"

rule test_genrich:
    input:
        "ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/genrich/G1/H2AZG1_mm.narrowPeaks"

include_prefix = "./rules/"
include:
    include_prefix + "run_alignment.smk"

