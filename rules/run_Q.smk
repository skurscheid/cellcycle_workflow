__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2019-04-09"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Rules for performing Q peak calling on ChIP-Seq data.

For use, include in your workflow.
"""

import os
import fnmatch
import snakemake
from snakemake.exceptions import MissingInputException

rule perform_q_peak_calling:
    version: 1
    params:
        cutoff = 0.005 # -log10(p-value)
    threads: 
        8
    input:
        chip = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{library}-{rep}.bam",
        input = "{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/INPUT{cycle}_{cycle}.bam"
    output:
        outdir = directory("{assayType}/{project}/{runID}/Q/callpeak/{reference_version}/{cycle}/{library}-{rep}/")
    shell:
    """
        if [ ! -d {output.outdir} ]; then mkdir -p {output.outdir}; fi &&\
        /home/skurscheid/Bioinformatics/Q/bin/Q --treatment-sample {input.chip}\
                                                --control-sample {input.input}\
                                                --out-prefix {output.outdir}{wildcards.library}-{wildcards.rep}\
                                                --thread-num {threads}\
                                                --p-value-cutoff {params.cutoff}
    """

rule run_q_all:
    input:
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/Q/callpeak/GRCh38_ensembl84/G1/{chip_library}-{rep}",
               chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["G1"]["ChIP"].keys()],
               rep = [1, 2]),
        expand("ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/Q/callpeak/GRCh38_ensembl84/M/{chip_library}-{rep}",
               chip_library = [x for x in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"]["M"]["ChIP"].keys()],
               rep = [1, 2])

rule run_q_single:
    input:
        "ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/Q/callpeak/GRCh38_ensembl84/G1/H2AZG1-1"