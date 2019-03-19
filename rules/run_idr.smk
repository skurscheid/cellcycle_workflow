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

rule run_idr:
    version: 
        "1"
    conda:
        "../envs/idr.yaml"
    params:
    input:
        rep1 = "{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{chip_library}-1/{chip_library}-1_peaks.narrowPeak",
        rep2 = "{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{chip_library}-2/{chip_library}-2_peaks.narrowPeak"
    output: 
        results = "{assayType}/{project}/{runID}/idr/pairwise/{reference_version}/{cycle}/{chip_library}_idr.narrowPeak"
    shell:
        """
            idr --samples {input.rep1} {input.rep2}\
                --input-file-type narrowPeak\
                --output-file {output.results}\
                --output-file-type narrowPeak\
                --log-output-file {output.results}/log.txt\
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
        idr_peaks = "{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{chip_library}_idr.bed",
        other_peaks = "{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{chip_library}_other.bed"
    script:
        "../scripts/extract_idr_peaks.py"

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
        bigwig_file = "{assayType}/{project}/{runID}/deepTools/bamCompare/{reference_version}/{cycle}/{chip_library}-{rep}_RPKM.bw"
    output:
        matrix = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{reference_version}/{cycle}/{chip_library}-{rep}.gz",
        bed = "{assayType}/{project}/{runID}/deepTools/computeMatrix/{reference_version}/{cycle}/{chip_library}-{rep}.bed",
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
        averageTypeSummaryPlot = "std"
    input:
        matrix = rules.compute_peaks_matrix_per_sample.output.matrix
    output:
        pdf = "{assayType}/{project}/{runID}/deepTools/plotHeatmap/{reference_version}/{cycle}/{chip_library}-{rep}.pdf"
    shell:
        """
            plotHeatmap --matrixFile {input.matrix}\
                        --outFileName {output.pdf}\
                        --averageTypeSummaryPlot {params.averageTypeSummaryPlot}\
                        --refPointLabel "Peak center"\
                        --perGroup
        """

rule upload_plots:
    version:
         "1"
    input:
        rules.plot_peaks_per_sample.output.pdf
    output:
        AS.remote("experiment/{assayType}/{project}/{runID}/deepTools/plotHeatmap/{reference_version}/{cycle}/{chip_library}-{rep}.pdf")
    run:
        shell("mv {input} {output}")
    
