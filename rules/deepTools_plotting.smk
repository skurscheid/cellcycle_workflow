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
        fn.append("/".join((path, "_".join((i, wildcards["operation"] + ".bw")))))
    return(fn)


def cli_parameters_normalization(wildcards):
    if wildcards["norm"] == "RPKM":
        a = "--normalizeUsingRPKM"
    elif wildcards["norm"] == "1xcoverage":
        a = " ".join(("--normalizeTo1x", config["references"][REF_GENOME]["effectiveSize"]))
    return(a)

def cli_parameters_computeMatrix(wildcards):
    params = config["program_parameters"]["deepTools"]["computeMatrix"][wildcards["subcommand"]][wildcards["region_type"]]
    return(params)

rule computeMatrix_scaled:
    conda:
        "../envs/deeptools.yaml"
    version:
        "1"
    params:
        cli_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        32
    input:
        file = get_computeMatrix_input,
        region = lambda wildcards: config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{subcommand}/{reference_version}/{region_type}/{region}/matrix_{operation}.gz"
    shell:
        """
        computeMatrix scale-regions --numberOfProcessors {threads} \
                                    --smartLabels \
                                    --missingDataAsZero \
                                    {params.cli_parameters}\
                                    --regionsFileName {input.region} \
                                    --scoreFileName {input.file} \
                                    --outFileName {output.matrix_gz}
        """

rule computeMatrix_scaled_input_normalised:
    conda:
        "../envs/deeptools.yaml"
    version:
        "1"
    params:
        cli_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        32
    input:
        M = lambda wildcards: expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bamCoverage/log2/GRCh38_ensembl84/M/{library}.bw",\
                                        library = list(config["samples"][wildcards["assayType"]]["conditions"][wildcards["runID"]]["M"]["ChIP"].values())),
        G1 = lambda wildcards: expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bamCoverage/log2/GRCh38_ensembl84/G1/{library}.bw",\
                                        library = list(config["samples"][wildcards["assayType"]]["conditions"][wildcards["runID"]]["G1"]["ChIP"].values())),
        region = lambda wildcards: config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards["region"]]
    output:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{subcommand}/{reference_version}/{region_type}/{region}/matrix_{operation}.gz"
    shell:
        """
        computeMatrix scale-regions --numberOfProcessors {threads} \
                                    --smartLabels \
                                    --missingDataAsZero \
                                    {params.cli_parameters}\
                                    --regionsFileName {input.region} \
                                    --scoreFileName {input.M} {input.G1) \
                                    --outFileName {output.matrix_gz}
        """

rule computeMatrix_refPoint:
    conda:
        "../envs/deeptools.yaml"
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
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/reference-point/{reference_version}/{region_type}/{region}/matrix_{operation}.gz"
    shell:
        """
            computeMatrix reference-point --numberOfProcessors {threads} \
                                                             --smartLabels \
                                                             --missingDataAsZero \
                                                             --upstream {params.upstream} \
                                                             --downstream {params.downstream} \
                                                             --regionsFileName {input.region} \
                                                             --scoreFileName {input.file} \
                                                             --outFileName {output.matrix_gz}
        """


rule plotProfile:
    conda:
        "../envs/deeptools.yaml"
    version:
        "1"
    params:
        dpi = 300,
        averageType = "mean",
        plotType = "se",
        plotTitle = "\"Mean coverage, scaled\"",
        numPlotsPerRow = 4
    threads:
        1
    input:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{subcommand}/{reference_version}/{region_type}/{region}/matrix_{operation}.gz",
    output:
        pdf =  "{assayType}/{project}/{runID}/deepTools/plotProfile/{subcommand}/{reference_version}/{region_type}/{region}/matrix_{operation}.pdf"
    shell:
        """
            plotProfile --matrixFile {input.matrix_gz}\
                        --outFileName {output.pdf}\
                        --dpi {params.dpi}\
                        --averageType {params.averageType}\
                        --plotType {params.plotType}\
                        --plotTitle {params.plotTitle}\
                        --numPlotsPerRow {params.numPlotsPerRow}
        """

rule plotHeatmap:
    conda:
        "../envs/deeptools.yaml"
    version:
        "1"
    params:
        dpi = 300,
        averageTypeSummaryPlot = "mean",
        plotTitle = "\"Mean coverage, scaled\"",
        numPlotsPerRow = 4
    threads:
        1
    input:
        matrix_gz = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{subcommand}/{reference_version}/{region_type}/{region}/matrix_{operation}.gz",
    output:
        pdf =  "{assayType}/{project}/{runID}/deepTools/plotHeatmap/{subcommand}/{reference_version}/{region_type}/{region}/matrix_{operation}.pdf"
    shell:
        """
            plotHeatmap --matrixFile {input.matrix_gz}\
                        --outFileName {output.pdf}\
                        --dpi {params.dpi}\
                        --averageTypeSummaryPlot {params.averageTypeSummaryPlot}\
                        --plotTitle {params.plotTitle}
        """


rule testPlot:
    input:
        "ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/scale-regions/GRCh38_ensembl84/repeats/LINEs_sample_50k/matrix_log2.pdf"
