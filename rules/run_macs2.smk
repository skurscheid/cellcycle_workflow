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
from snakemake.remote.AzureStorage import RemoteProvider as AzureRemoteProvider

# setup Azure Storage for remote access
AS = AzureRemoteProvider(account_name= config["account_name"], 
    account_key=config["account_key"])

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
        genomeSize = "mm",
        name = lambda wildcards: wildcards.treatment
    input:
        treatment = AS.remote("experiment/{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{treatment}-{rep}.bam"),
        control = AS.remote("experiment/{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam")
    output:
        outDir = directory("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-{rep}")
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
