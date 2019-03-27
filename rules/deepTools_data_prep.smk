__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2019-02-25"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
from snakemake.io import expand
import os

"""
Rules for running deepTools analysis on ChIP-Seq data
For usage, include this in your workflow.
"""

def get_computeMatrix_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayType"],
                     wildcards["project"],
                     wildcards["runID"],
                     "deepTools",
                     "bamCoverage",
                     wildcards["reference_version"]))
    for i in config["samples"][wildcards["assayType"]][wildcards["project"]][wildcards["runID"]].keys():
        fn.append("/".join((path, "_".join((i, wildcards["suffix"] + ".bw")))))
    return(fn)

def cli_parameters_normalization(wildcards):
    if wildcards["norm"] == "RPKM":
        a = "--normalizeUsingRPKM"
    elif wildcards["norm"] == "1xcoverage":
        a = " ".join(("--normalizeTo1x", config["references"][REF_GENOME]["effectiveSize"]))
    return(a)

def cli_parameters_bamCoverage(wildcards):
    a = config["program_parameters"]["deepTools"]["bamCoverage"]["normal"]
    b = str()
    for (key, val) in a.items():
        if val == " ":
            f = key + " "
            b = b + f
        else:
            f = key + "=" + val + " "
            b = b + f
    return(b.rstrip())

def get_input_library(wildcards):
    libraries = config["samples"][wildcards["assayType"]]["conditions"][wildcards["runID"]][wildcards["condition"]]["Input"].keys()
    b = expand("{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{library}_{condition}.bam",
               assayType = wildcards["assayType"],
               project = wildcards["project"],
               runID = wildcards["runID"],
               reference_version = wildcards["reference_version"],
               library = libraries,
               condition = wildcards["condition"])
    return(b)

rule bamCoverage_normal:
    version:
        1
    conda:
        "envs/deeptools.yaml"
    params:
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        outFileFormat = "bigwig",
        binSize = 10,
        smoothLength = 30,
        normalizeUsing = "RPKM",
    threads:
        32
    input:
        bam = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{library}.bam",
        index = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{library}.bam.bai"
    output:
        bigwig = "{assayType}/{project}/{runID}/deepTools/bamCoverage/{reference_version}/{cycle}/{library}_RPKM.bw"
    shell:
        """
        bamCoverage --bam {input.bam} \
                    --outFileName {output.bigwig} \
                    --outFileFormat {params.outFileFormat} \
                    --binSize {params.binSize} \
                    --smoothLength {params.smoothLength}\
                    --numberOfProcessors {threads} \
                    --normalizeUsing RPKM \
                    --extendReads \
                    --ignoreForNormalization {params.ignore}
        """

rule bamCompare:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    params:
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        outFileFormat = "bigwig",
        binSize = 25,
        smoothLength = 50,
        normalizeUsing = "RPKM",
    threads:
        8
    input:
        chip_library_bam = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{chip_library}-{rep}.bam",
        chip_library_index = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{chip_library}-{rep}.bam.bai",
        control_bam = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.bam",
        control_index = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.bam.bai"
    output:
        bigwig = "{assayType}/{project}/{runID}/deepTools/bamCompare/{reference_version}/{cycle}/{chip_library}-{rep}_RPKM.bw"
    shell:
        """
            bamCompare --bamfile1 {input.chip_library_bam}\
                       --bamfile2 {input.control_bam}\
                       --outFileName {output.bigwig}\
                       --outFileFormat bigwig\
                       --normalizeUsingRPKM\
                       --ignoreForNormalization {params.ignore}\
                       --smoothLength {params.smoothLength}\
                       --binSize {params.binSize} \
                       --numberOfProcessors {threads}
        """

rule subtractWT:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    params:
    threads:
        8
    input:
        chip_library_bw = "{assayType}/{project}/{runID}/deepTools/bamCoverage/{reference_version}/{chip_library}-{rep}_RPKM.bw",
        wt_bw = "{assayType}/{project}/{runID}/deepTools/bamCoverage/{reference_version}/WT{cycle}-{rep}_RPKM.bw",
    output:
        bigwig = "{assayType}/{project}/{runID}/deepTools/bigwigCompare/{operation}/{reference_version}/{cycle}/{chip_library}-{rep}.bw"
    shell:
        """
            bigwigCompare --bigwig1 {input.chip_library_bw}\
                       --bigwig2 {input.wt_bw}\
                       --outFileName {output.bigwig}\
                       --numberOfProcessors {threads}\
                       --operation {wildcards.operation}
        """

rule subtractOneSample:
    input:
      "ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bigwigCompare/subtract/GRCh38_ensembl84/G1/ACTR6G1-1.bw"

rule runSubtractWT:
    input:
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bigwigCompare/subtract/GRCh38_ensembl84/G1/{library}-{rep}.bw",
               library = ["ACTR6M", "ANP32M", "H2AM", "H2AZM", "TIP60M", "YL1M"],
               rep = ["1", "2"]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bigwigCompare/subtract/GRCh38_ensembl84/M/{library}-{rep}.bw",
               library = ["ACTR6G1", "ANP32G1", "H2AG1", "H2AZG1", "TIP60G1", "YL1G1"],
               rep = ["1", "2"])
