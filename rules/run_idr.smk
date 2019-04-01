__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2019-03-14"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for performing irreproducibility of discovery rate (IDR) analysis on peak calls.

For use, include in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException
from snakemake.remote.AzureStorage import RemoteProvider as AzureRemoteProvider

# setup Azure Storage for remote access
account_key=os.environ['AZURE_KEY']
account_name=os.environ['AZURE_ACCOUNT']
AS = AzureRemoteProvider(account_name=account_name, account_key=account_key)

rule run_idr:
    version: 
        "1"
    conda:
        "../envs/idr.yaml"
    params:
    input:
        rep1 = "{assayType}/{project}/{runID}/macs2/{command}/{reference_version}/{cycle}/{chip_library}-1/{chip_library}_peaks.narrowPeak",
        rep2 = "{assayType}/{project}/{runID}/macs2/{command}/{reference_version}/{cycle}/{chip_library}-2/{chip_library}_peaks.narrowPeak"
    output: 
        results = "{assayType}/{project}/{runID}/idr/{command}/{reference_version}/{cycle}/{chip_library}_idr.narrowPeak"
    shell:
        """
            idr --samples {input.rep1} {input.rep2}\
                --input-file-type narrowPeak\
                --output-file {output.results}\
                --output-file-type narrowPeak\
                --log-output-file {output.results}.log\
                --use-best-multisummit-IDR\
        """

rule extract_idr_peaks:
    version:
        "1"
    conda:
        "../envs/pandas.yaml"
    params:
        globalIDRCutoff = 1,
        signalValue = 5
    input:
        idr_file = rules.run_idr.output.results
    output:
        idr_peaks = "{assayType}/{project}/{runID}/idr/{command}/BEDs/{reference_version}/{cycle}/{chip_library}_idr.bed",
        other_peaks = "{assayType}/{project}/{runID}/idr/{command}/BEDs/{reference_version}/{cycle}/{chip_library}_other.bed"
    script:
        "../scripts/extract_idr_peaks.py"

rule bigWigCompare_vs_Input:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    params:
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        outFileFormat = "bigwig",
        binSize = 25,
        smoothLength = 50,
    threads:
        8
    input:
        chip_bw = "{assayType}/{project}/{runID}/deepTools/bigwigCompare/subtract/{reference_version}/{cycle}/{chip_library}{cycle}-{rep}_RPKM.bw",
        control_bw = "{assayType}/{project}/{runID}/deepTools/bamCoverage/{reference_version}/{cycle}/INPUT{cycle}_RPKM.bw"
    output:
        bigwig = "{assayType}/{project}/{runID}/deepTools/bamCompare/{reference_version}/{cycle}/{chip_library}-{rep}_log2.bw"
    shell:
        """
            bamCompare --bamfile1 {input.chip_library_bam}\
                       --bamfile2 {input.control_bam}\
                       --outFileName {output.bigwig}\
                       --outFileFormat bigwig\
                       --ignoreForNormalization {params.ignore}\
                       --smoothLength {params.smoothLength}\
                       --binSize {params.binSize} \
                       --numberOfProcessors {threads}
        """

rule compute_peaks_matrix_per_sample:
    version:
        "1"
    conda:
        "../envs/deeptools.yaml"
    params:
        beforeRegionStartLength = 1000,
        afterRegionStartLength = 1000
    threads:
        16
    input:
        idr_peaks = rules.extract_idr_peaks.output.idr_peaks,
        other_peaks = rules.extract_idr_peaks.output.other_peaks,
        bigwig_file = "{assayType}/{project}/{runID}/deepTools/bamCompare/{reference_version}/{cycle}/{chip_library}-{rep}_log2.bw"
    output:
        matrix = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{command}/{reference_version}/{cycle}/{chip_library}-{rep}.gz",
        bed = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{command}/{reference_version}/{cycle}/{chip_library}-{rep}.bed",
    shell:
        """
            computeMatrix reference-point --scoreFileName {input.bigwig_file}\
                                          --regionsFileName {input.idr_peaks} {input.other_peaks}\
                                          --outFileName {output.matrix}\
                                          --outFileSortedRegions {output.bed}\
                                          --referencePoint center\
                                          --beforeRegionStartLength {params.beforeRegionStartLength}\
                                          --afterRegionStartLength {params.afterRegionStartLength}\
                                          --numberOfProcessors {threads}
        """
    
rule plot_peaks_per_sample:
    version:
        "1"
    conda:
        "../envs/deeptools.yaml"
    params:
        averageTypeSummaryPlot = "mean"
    input:
        matrix = rules.compute_peaks_matrix_per_sample.output.matrix
    output:
        pdf = "{assayType}/{project}/{runID}/deepTools/plotHeatmap/{command}/{reference_version}/{cycle}/{chip_library}-{rep}-{summaryPlotType}.pdf"
    shell:
        """
            plotHeatmap --matrixFile {input.matrix}\
                        --outFileName {output.pdf}\
                        --averageTypeSummaryPlot {wildcards.summaryPlotType}\
                        --refPointLabel "Peak center"\
                        --perGroup
        """

rule allplots:
    input:
         expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/{command}/GRCh38_ensembl84/G1/{chip_library}-{rep}-{summaryPlotType}.pdf",
                command = "callpeak_combined_controls",
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
                rep = ["1", "2"],
                summaryPlotType = "mean"),
         expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/{command}/GRCh38_ensembl84/M/{chip_library}-{rep}-{summaryPlotType}.pdf",
                command = "callpeak_combined_controls",
                chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
                rep = ["1", "2"],
                summaryPlotType = "mean")
