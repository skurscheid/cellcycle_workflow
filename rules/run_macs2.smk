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
        treatment = AS.remote("experiment/{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{treatment}-{rep}.{ext}")
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{treatment}-{rep}.{ext}"
    run:
        shell("cp {input.treatment} {output}")

rule download_control:
    version:
        "1"
    input:
        control = AS.remote("experiment/{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.{ext}")
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.{ext}"
    run:
        shell("cp {input.control} {output}")

rule index_control:
    version:
        "1"
    conda:
        "../envs/samtools.yaml"
    input:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.bam"
    output:
        "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.bam.bai"
    shell:
        "samtools index {input} {output}"

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
        treatment_bam = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{treatment}-{rep}.bam",
        treatment_index = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/{treatment}-{rep}.bam.bai",
        control_bam = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.bam",
        control_index = "{assayType}/{project}/{runID}/transfer/down/{reference_version}/{cycle}/INPUT{cycle}_{cycle}.bam.bai"
    output:
        outDir = directory("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{treatment}-{rep}")
    shell:
        """
            macs2 callpeak -f {params.fileType} \
                           -g hs\
                           -t {input.treatment_bam}\
                           -c {input.control_bam}\
                           -n {params.name}\
                           --outdir {output.outDir}\
                           --call-summits\
                           --qvalue {params.qvalCutoff}\
                           --bdg\
                           --trackline\
                           --tempdir /data/tmp
        """
