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
home = os.environ['HOME']
REF_GENOME = config["references"]["active"]
REF_VERSION = config["references"][REF_GENOME]["version"][0]


rule bowtie2_pe:
    conda:
        "../envs/bam.yaml"
    version:
        "1"
    params:
        max_in = config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index = config["references"][REF_GENOME]["bowtie2"][REF_VERSION]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    input:
        trimmed_read1 = "{assayType}/{project}/{runID}/fastp/trimmed/{library}.end1.fastq.gz",
        trimmed_read2 = "{assayType}/{project}/{runID}/fastp/trimmed/{library}.end2.fastq.gz"
    output:
        temp("{assayType}/{project}/{runID}/bowtie2/{library}_pe.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {threads}\
            --rg-id '{wildcards.library}' \
            --rg 'LB:{wildcards.library}' \
            --rg 'SM:{wildcards.library}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.trimmed_read1} \
            -2 {input.trimmed_read2} \
            | samtools view -Sb - > {output}
        """

rule bowtie2_pe_multimap:
    conda:
        "../envs/bam.yaml"
    version:
        "1"
    params:
        max_in = config["program_parameters"]["bt2_params"]["max_insert"],
        bt2_index = config["references"][REF_GENOME]["bowtie2"][REF_VERSION]
    threads:
        lambda wildcards: int(str(config["program_parameters"]["bt2_params"]["threads"]).strip("['']"))
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_mm_bowtie2_pe_multimap_times.tsv"
    input:
        trimmed_read1 = "{assayType}/{project}/{runID}/fastp/trimmed/{library}.end1.fastq.gz",
        trimmed_read2 = "{assayType}/{project}/{runID}/fastp/trimmed/{library}.end2.fastq.gz"
    output:
        temp("{assayType}/{project}/{runID}/bowtie2/{library}_mm.bam")
    shell:
        """
            bowtie2 \
            -x {params.bt2_index}\
            --no-mixed \
            --no-discordant \
            --maxins {params.max_in} \
            --threads {threads}\
            -k 5\
            --rg-id '{wildcards.library}' \
            --rg 'LB:{wildcards.library}' \
            --rg 'SM:{wildcards.library}' \
            --rg 'PL:Illumina' \
            --rg 'PU:NA' \
            -1 {input.trimmed_read1} \
            -2 {input.trimmed_read2} \
            | samtools view -Sb - > {output}
        """


rule bam_stats:
    conda:
        "../envs/bam.yaml"
    version:
        "1"
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_stats_times.tsv"
    input:
        "{assayType}/{project}/{runID}/bowtie2/{library}_{suffix}.bam"
    output:
        "{assayType}/{project}/{runID}/samtools/flagstat/{library}_{suffix}.bam.stats.txt"
    shell:
        """
            samtools flagstat {input} > {output}
        """


# rules
rule bam_quality_filter:
    conda:
        "../envs/bam.yaml"
    version:
        "1.0"
    params:
        qual = config["alignment_quality"]
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_quality_filter_times.tsv"
    input:
        "{assayType}/{project}/{runID}/bowtie2/{library}_{suffix}.bam"
    output:
        temp("{assayType}/{project}/{runID}/samtools/quality_filtered/{library}_{suffix}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output}"


rule bam_sort:
    conda:
        "../envs/bam.yaml"
    version:
        "1.0"
    threads:
        4
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_sort_times.tsv"
    input:
        rules.bam_quality_filter.output
    output:
        temp("{assayType}/{project}/{runID}/samtools/sort/{library}_{suffix}.bam")
    shell:
        "samtools sort -@ {threads} {input} -T {wildcards.library}.sorted -o {output}"

rule bam_mark_duplicates:
    conda:
        "../envs/bam.yaml"
    version:
        "1.0"
    params:
        qual = config["alignment_quality"],
        picard = home + config["program_parameters"]["picard_tools"]["jar"],
        temp = home + config["temp_dir"]
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_mark_duplicates_times.tsv"
    input:
        rules.bam_sort.output
    output:
        out= "{assayType}/{project}/{runID}/picardTools/MarkDuplicates/{library}_{suffix}.bam",
        metrics = "{assayType}/{project}/{runID}/picardTools/MarkDuplicates/{library}_{suffix}.metrics.txt"
    shell:
        """
            picard MarkDuplicates\
            INPUT={input}\
            OUTPUT={output.out}\
            ASSUME_SORTED=TRUE\
            REMOVE_DUPLICATES=TRUE\
            METRICS_FILE={output.metrics}
        """

rule bam_sortn:
    conda:
        "../envs/bam.yaml"
    version:
        "1.0"
    threads:
        4
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_sortn_times.tsv"
    input:
        rules.bam_mark_duplicates.output
    output:
        "{assayType}/{project}/{runID}/samtools/sortn/{library}_{suffix}.bam"
    shell:
        "samtools sort -n -@ {threads} {input} -T {wildcards.library}.sorted -o {output}"

rule bam_index:
    conda:
        "../envs/bam.yaml"
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_index_times.tsv"
    input:
        rules.bam_mark_duplicates.output
    output:
        "{assayType}/{project}/{runID}/picardTools/MarkDuplicates/{library}_{suffix}.bam.bai"
    shell:
        "samtools index {input} {output}"

rule bam_insert_size:
    conda:
        "../envs/bam.yaml"
    version:
        "1.0"
    params:
        picard = home + config["program_parameters"]["picard_tools"]["jar"],
        temp = home + config["temp_dir"]
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/{library}_{suffix}_bam_insert_size_times.tsv"
    input:
        rules.bam_mark_duplicates.output
    output:
        metrics = "{assayType}/{project}/{runID}/picardTools/CollectInsertSizeMetrics/{library}_{suffix}.insert_size_metrics.txt",
        histogram = "{assayType}/{project}/{runID}/picardTools/CollectInsertSizeMetrics/{library}_{suffix}.histogram.pdf"
    shell:
        """
            picard CollectInsertSizeMetrics \
            INPUT={input}\
            OUTPUT={output.metrics}\
            H={output.histogram}
        """
