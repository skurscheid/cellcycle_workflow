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
        rep1 = "{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-1/{treatment}_peaks.narrowPeak",
        rep2 = "{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-2/{treatment}_peaks.narrowPeak"
    output:
        results = "{assayType}/{project}/{runID}/idr/pairwise/{reference_version}/{cycle}/{treatment}/{treatment}_idr.narrowPeak",
        log = "{assayType}/{project}/{runID}/idr/pairwise/{reference_version}/{cycle}/{treatment}/{treatment}_idr.log"
    shell:
        """
            idr --samples {input.rep1} {input.rep2}\
                --input-file-type narrowPeak\
                --output-file {output.results}\
                --output-file-type narrowPeak\
                --log-output-file {output.log}\
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
        idr_peaks = "{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{treatment}/{treatment}_idr.bed",
        other_peaks = "{assayType}/{project}/{runID}/idr/BEDs/{reference_version}/{cycle}/{treatment}/{treatment}_other.bed"
    script:
        "../scripts/extract_idr_peaks.py"

    
