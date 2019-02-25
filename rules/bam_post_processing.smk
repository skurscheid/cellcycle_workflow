# bam post-processing

from snakemake.exceptions import MissingInputException
import os
import pandas as pd
home = os.environ['HOME']

"""
Rules for post-processing BAM files
For usage, include this in your workflow.
"""

def getSampleInformation(wildcards):
    tab = pd.read_csv(config["samples"][wildcards["assayType"]]["libraries"][wildcards["project"]][wildcards["runID"]]["input_libraries"]).set_index("condition")
    return(tab)

def getReplicates(wildcards):
    tab = getSampleInformation(wildcards)
    fn =[]
    fn = expand("{assayType}/{project}/{runID}/bamtools/merge/{reference_version}/{input}.bam",
                assayType = wildcards["assayType"],
                project = wildcards["project"],
                runID = wildcards["runID"],
                tool = wildcards["tool"],
                command = wildcards["command"],
                reference_version = wildcards["reference_version"],
                input = input_libraries.loc[wildcards["condition"], "sample"])
    return(fn)

rule merge_input:
    version:
        1
    conda:
        "envs/bamtools.yaml"
    params:
        conditions = config["samples"]["conditions"] 
    threads:
        8
    input:
        bamFiles = getReplicates
    output:
        merged = "{assayType}/{project}/{runID}/bamtools/merge/{reference_version}/{mergedLibrary}_{condition}.bam"
    shell:
        """
            bamtools -in {input.bamFiles}\
                     -out {output.merged}
        """


