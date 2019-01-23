__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-10"

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
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
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
        {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                           --outFileName {output.bigwig} \
                                           --outFileFormat {params.outFileFormat} \
                                           --binSize {params.binSize} \
                                           --smoothLength {params.smoothLength}\
                                           --numberOfProcessors {threads} \
                                           --normalizeUsing RPKM \
                                           --extendReads \
                                           --ignoreForNormalization {params.ignore}
        """

rule computeMatrix_scaled:
    version:
        "1"
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        regionBodyLength = 5000,
        beforeRegionStartLength = 2000,
        afterRegionStartLength = 2000,
        unscaled5prime = 350,
        unscaled3prime = 350
    threads:
        32
    input:
        file = get_computeMatrix_input,
        region = lambda wildcards: config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/scale-region/{reference_version}/{region}/matrix_{suffix}.gz"
    shell:
        """
        {params.deepTools_dir}/computeMatrix scale-regions --numberOfProcessors {threads} \
                                                           --smartLabels \
                                                           --missingDataAsZero \
                                                           --regionBodyLength {params.regionBodyLength} \
                                                           --beforeRegionStartLength {params.beforeRegionStartLength} \
                                                           --afterRegionStartLength {params.afterRegionStartLength} \
                                                           --unscaled5prime {params.unscaled5prime} \
                                                           --unscaled3prime {params.unscaled3prime} \
                                                           --regionsFileName {input.region} \
                                                           --scoreFileName {input.file} \
                                                           --outFileName {output.matrix_gz}
        """


rule computeMatrix_refPoint:
    version:
        "1"
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        upstream = 500,
        downstream = 500
    threads:
        32
    input:
        file = get_computeMatrix_input,
        region = lambda wildcards: config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/reference-point/{reference_version}/{region}/matrix_{suffix}.gz"
    shell:
        """
        {params.deepTools_dir}/computeMatrix reference-point --numberOfProcessors {threads} \
                                                             --smartLabels \
                                                             --missingDataAsZero \
                                                             --upstream {params.upstream} \
                                                             --downstream {params.downstream} \
                                                             --regionsFileName {input.region} \
                                                             --scoreFileName {input.file} \
                                                             --outFileName {output.matrix_gz}
        """


rule plotProfile:
    version:
        "1"
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        dpi = 300,
        averageType = "mean",
        plotType = "se",
        plotTitle = "\"Mean coverage, all genes, scaled\"",
        numPlotsPerRow = 4
    threads:
        1
    input:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{subcommand}/{reference_version}/{region}/matrix_{suffix}.gz",
    output:
        pdf =  "{assayType}/{project}/{runID}/deepTools/plotProfile/{subcommand}/{reference_version}/{region}_{suffix}.pdf"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz}\
                                               --outFileName {output.pdf}\
                                               --dpi {params.dpi}\
                                               --averageType {params.averageType}\
                                               --plotType {params.plotType}\
                                               --plotTitle {params.plotTitle}\
                                               --numPlotsPerRow {params.numPlotsPerRow}
        """
