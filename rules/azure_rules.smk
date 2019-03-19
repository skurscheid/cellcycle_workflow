__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for performing macs2 peak calling on ChIP-Seq data.

For use, include in your workflow.
"""

import os
import fnmatch
import snakemake
from snakemake.exceptions import MissingInputException
from snakemake.remote.AzureStorage import RemoteProvider as AzureRemoteProvider

# setup Azure Storage for remote access
account_key=os.environ['AZURE_KEY']
account_name=os.environ['AZURE_ACCOUNT']
AS = AzureRemoteProvider(account_name=account_name, account_key=account_key)

# set local variables
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]

rule download_treatment:
    version:
        "1"
    input:
        treatment = AS.remote("experiment/{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{chip_library}-{rep}.bam")
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{chip_library}-{rep}.bam"
    run:
        shell("mv {input.treatment} {output}")

rule index_treatment:
    version:
        "1"
    conda:
        "../envs/samtools.yaml"
    input:
        rules.download_treatment.output
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{chip_library}-{rep}.bam.bai"
    shell:
        "samtools index {input} {output}"

rule download_control:
    version:
        "1"
    input:
        control = AS.remote("experiment/{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam")
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/INPUT{cycle}_{cycle}.bam"
    run:
        shell("mv {input.control} {output}")

rule index_control:
    version:
        "1"
    conda:
        "../envs/samtools.yaml"
    input:
        rules.download_control.output
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/INPUT{cycle}_{cycle}.bam.bai"
    shell:
        "samtools index {input} {output}"