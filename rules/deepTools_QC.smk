__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-08-07"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools QC on ChIP-Seq data
For usage, include this in your workflow.
"""

# functions
def librariesPerCondition(wildcards):
    libs = []
    libs = expand("{assayType}/{project}/{runID}/samtools/rmdup/{reference_version}/{library}{replicate}.{suffix}",
                  assayType = "ChIP-Seq",
                  project = PROJECT_ID,
                  runID = RUN_ID,
                  reference_version = REF_VERSION,
                  library = [z \
                              for y in config["samples"]["ChIP-Seq"]["conditions"][RUN_ID][wildcards["condition"]].keys() \
                                  for z in config["samples"]["ChIP-Seq"]["conditions"][RUN_ID][wildcards["condition"]][y]],
                  replicate = ["-1", "-2"],
                  suffix = ["bam"])
    return(libs)

def sampleLabelsPerCondition(wildcards):
    labels = []
    labels = expand("{library}{replicate}",
                    library = [z \
                                for y in config["samples"]["ChIP-Seq"]["conditions"][RUN_ID][wildcards["condition"]].keys() \
                                    for z in config["samples"]["ChIP-Seq"]["conditions"][RUN_ID][wildcards["condition"]][y]],
                    replicate = ["-1", "-2"])
    return(labels)

rule multiBamSummary:
    version:
        "2"
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        sampleLabels = sampleLabelsPerCondition
    threads:
        32
    input:
        librariesPerCondition
    output:
        npz = "{assayType}/{project}/{runID}/deepTools/multiBamSummary/{reference_version}/{condition}/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input} \
                                                        --numberOfProcessors {threads} \
                                                        --labels {params.sampleLabels} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """


rule plotCorrelation_heatmap:
    version:
        "2"
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = "Correlation heatmap - read counts",
        sampleLabels = sampleLabelsPerCondition
    input:
        npz = rules.multiBamSummary.output.npz
    output:
        png = "{assayType}/{project}/{runID}/deepTools/plotCorrelation/{reference_version}/{condition}/heatmap_SpearmanCorr_readCounts.png",
        tab = "{assayType}/{project}/{runID}/deepTools/plotCorrelation/{reference_version}/{condition}/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
            {params.deepTools_dir}/plotCorrelation --corData {input.npz} \
                                                   --corMethod spearman \
                                                   --skipZeros \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.sampleLabels} \
                                                   --whatToPlot heatmap \
                                                   --colorMap RdYlBu \
                                                   --plotNumbers \
                                                   -o {output.png} \
                                                   --outFileCorMatrix {output.tab}
        """

rule plotPCA:
    version:
        "2"
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = "PCA - read counts",
        sampleLabels = sampleLabelsPerCondition
    input:
        npz = rules.multiBamSummary.output.npz
    output:
        png = "{assayType}/{project}/{runID}/deepTools/plotPCA/{reference_version}/{condition}/PCA_readCounts.png"
    shell:
        """
            {params.deepTools_dir}/plotPCA --corData {input.npz} \
                                           --labels {params.sampleLabels} \
                                           --plotFile {output.png} \
                                           --plotTitle "{params.plotTitle}"
        """

rule bamPEFragmentSize:
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = "BAM PE fragment size",
        sampleLabels = sampleLabelsPerCondition
    threads:
        32
    input:
        librariesPerCondition
    output:
        "{assayType}/{project}/{runID}/deepTools/bamPEFragmentSize/{reference_version}/{condition}/histogram.png"
    shell:
        """
            {params.deepTools_dir}/bamPEFragmentSize --bamfiles {input} \
                                                     --numberOfProcessors {threads} \
                                                     --samplesLabel {params.sampleLabels}\
                                                     --histogram {output}
        """


rule plotFingerprint:
    params:
        deepTools_dir = home + config["program_parameters"]["deepTools"]["deepTools_dir"],
        plotTitle = "BAM PE fingerprint",
        sampleLabels = sampleLabelsPerCondition
    threads:
        32
    input:
        librariesPerCondition
    output:
        "{assayType}/{project}/{runID}/deepTools/plotFingerprint/{reference_version}/{condition}/fingerprints.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input} \
                                                   --numberOfProcessors {threads} \
                                                   --centerReads \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.sampleLabels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """
