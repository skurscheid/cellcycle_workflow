__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2019-02-25"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os
home = os.environ['HOME']

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
        bam = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{library}.bam",
        index = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{library}.bam.bai"
    output:
        bigwig = "{assayType}/{project}/{runID}/deepTools/bamCoverage/{reference_version}/{library}_RPKM.bw"
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
        "envs/deeptools.yaml"
    params:
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        outFileFormat = "bigwig",
        binSize = 25,
        smoothLength = 50,
        normalizeUsing = "RPKM",
    threads:
        8
    input:
        chip = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{chip}.bam"
        input = "{assayType}/{project}/{runID}/bamtools/merge/{reference_version}/{mergedLibrary}_{condition}.bam"
    output:
        bigwig = "{assayType}/{project}/{runID}/deepTools/bamCompare/{reference_version}/{chip_library}_{condition}_RPKM.bw"
    shell:
        """
            bamCompare --bamfile1 {input.chip}\
                       --bamfile2 {input.input}\
                       --outFileName {output.bigwig}\
                       --outFileFormat bigwig\
                       --normalizeUsing {params.normalizeUsing}\
                       --ignoreForNormalization {params.ignore}\
                       --smoothLength {params.smoothLength}\
                       --binSize {params.binSize}   
        """


