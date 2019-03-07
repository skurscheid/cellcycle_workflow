__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for aligning paired-end reads using bowtie2.

For use, include in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set local variables
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]

rule macs2_callpeak:
    conda:
        "../envs/macs2.yaml"
    version:
        "1"
    params:
        seed = "1234",
        fileType = "BAMPE",
        qvalCutoff = 0.99,
        genomeSize = "mm"
        name = lambda wildcards: wildcards.treatment
    input:
        treatment = "{assayType}/{project}/{runID}/samtools/rmdup/{treatment}.bam",
        control = "{assayType}/{project}/{runID}/samtools/rmdup/Input.bam"
    output:
        outDir = "{assayType}/{project}/{runID}/macs2/callpeak/{treatment}",
        bed = "{assayType}/{project}/{runID}/macs2/callpeak/{treatment}_summits.bed",
        xls = "{assayType}/{project}/{runID}/macs2/callpeak/{treatment}_peaks.xls"
    shell:
        """
            macs2 callpeak -f {params.fileType} \
                           -g hs\
                           -t {input.treatment}\
                           -c {input.control}
                           -n {params.name}\
                           --outdir {output.outDir}\
                           --call-summits\
                           ---qvalue {params.qvalCutoff}\
                           --bdg\
                           --trackline
        """
