rule genrich:
    conda:
        "../envs/bam.yaml"
    version:
        "1.0"
    threads:
        1
    params:
        genrich_path = "/short/kv78/software/Genrich",
        Ns_file = "/short/kv78/References/Genomes/Homo_sapiens/GRCh38_ensembl84/primary/Homo_sapiens.GRCh38.dna.primary_assembly.Ns.bed",
        min_auc = 10,
        max_qval = 0.1,
        additional_cli = genrich_parameters
    benchmark:
        "{assayType}/{project}/{runID}/benchmarks/genrich/{cycle}/{sample}_{suffix}_times.tsv"
    log:
        "{assayType}/{project}/{runID}/logs/genrich/{cycle}/{sample}_{suffix}.txt"
    input:
        chip = chip_replicates,
        input = input_replicates
    output:
        "{assayType}/{project}/{runID}/genrich/{cycle}/{sample}_{suffix}.narrowPeaks"
    shell:
        """
            {params.genrich_path}/Genrich {params.additional_cli}\
                                          -E {params.Ns_file}\
                                          -a {params.min_auc}\
                                          -q {params.max_qval}\
                                          -t {input.chip[0]},{input.chip[1]}\
                                          -c {input.input[0]},{input.input[1]}\
                                          -o {output} 1>>{log} 2>>{log}
        """
