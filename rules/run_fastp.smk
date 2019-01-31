__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

import os
import fnmatch
from snakemake.exceptions import MissingInputException

home = os.environ['HOME']

singularity: "docker://skurscheid/snakemake_baseimage:0.1"

rule run_fastp:
    conda:
        "envs/fastp.yaml"
    version:
        "1.2"
    threads:
        4
    input:
        read1 = lambda wildcards: dataDir + wildcards["assayType"] + "/" + wildcards["project"] + "/" + wildcards["runID"] + "/raw_data/" + config["samples"][wildcards["assayType"]][wildcards["project"]][wildcards["runID"]][wildcards["library"]][0],
        read2 = lambda wildcards: dataDir + wildcards["assayType"] + "/" + wildcards["project"] + "/" + wildcards["runID"] + "/raw_data/" + config["samples"][wildcards["assayType"]][wildcards["project"]][wildcards["runID"]][wildcards["library"]][1]
    output:
        trimmed_read1 = "{assayType}/{project}/{runID}/fastp/trimmed/{library}.end1.fastq.gz",
        trimmed_read2 = "{assayType}/{project}/{runID}/fastp/trimmed/{library}.end2.fastq.gz",
        report_html = "{assayType}/{project}/{runID}/fastp/report/{library}.fastp.html",
        report_json = "{assayType}/{project}/{runID}/fastp/report/{library}.fastp.json"
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.trimmed_read1} -O {output.trimmed_read2} --html {output.report_html} --json {output.report_json} --thread {threads}"
