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

rule macs2_callpeak:
    conda:
        "../envs/macs2.yaml"
    version:
        "1"
    params:
        seed = "1234",
        fileType = "BAMPE",
        qvalCutoff = 0.99,
        genomeSize = "hs",
        name = lambda wildcards: wildcards.chip_library
    input:
        treatment_bam = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{chip_library}-{rep}.bam",
        treatment_index = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{chip_library}-{rep}.bam.bai",
        control_bam = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam",
        control_index = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam.bai"
    output:
        outDir = directory("{assayType}/{project}/{runID}/macs2/callpeak/{reference_version}/{cycle}/{chip_library}-{rep}")
    shell:
        """
            macs2 callpeak -f {params.fileType}\
                           -g {params.genomeSize}\
                           -t {input.treatment_bam}\
                           -c {input.control_bam}\
                           -n {params.name}\
                           --outdir {output.outDir}\
                           --call-summits\
                           --bdg\
                           --qvalue {params.qvalCutoff}\
                           --tempdir tmp
        """

rule macs2_callpeak_combined_controls:
    conda:
        "../envs/macs2.yaml"
    version:
        "1"
    params:
        seed = "1234",
        fileType = "BAMPE",
        qvalCutoff = 0.99,
        genomeSize = "hs",
        name = lambda wildcards: wildcards.chip_library
    input:
        treatment_bam = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{chip_library}-{rep}.bam",
        treatment_index = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{chip_library}-{rep}.bam.bai",
        control_bam = ["{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam", 
                       "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/WT{cycle}-1.bam",
                       "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/WT{cycle}-2.bam"],
        control_index = ["{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam.bai",
                         "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/WT{cycle}-1.bam.bai",
                         "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/WT{cycle}-2.bam.bai"]
    output:
        outDir = directory("{assayType}/{project}/{runID}/macs2/callpeak_combined_controls/{reference_version}/{cycle}/{chip_library}-{rep}")
    shell:
        """
            macs2 callpeak -f {params.fileType}\
                           -g {params.genomeSize}\
                           -t {input.treatment_bam}\
                           -c {input.control_bam}\
                           -n {params.name}\
                           --outdir {output.outDir}\
                           --call-summits\
                           --qvalue {params.qvalCutoff}\
                           --tempdir tmp
        """

rule run_macs2_callpeak_combined_controls:
    input:
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/macs2/callpeak_combined_controls/GRCh38_ensembl84/G1/{chip_library}-{rep}",
               chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
               rep = [1, 2]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/macs2/callpeak_combined_controls/GRCh38_ensembl84/M/{chip_library}-{rep}",
               chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
               rep = [1, 2])