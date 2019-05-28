#!/bin/bash
#PBS -P kv78
#PBS -l walltime=24:00:00
#PBS -l wd
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=1
#PBS -l mem=4G
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -N snakemake_master
#PBS -q express

source ~/.bashrc

/short/rl2/miniconda3/envs/snakemake/bin/snakemake -s /short/kv78/cellcycle_workflow/Snakefile test_genrich\
        --configfile /short/kv78/cellcycle_workflow/config.json\
        --use-conda\
        --jobs 32\
        -d /short/kv78/cellcycle\
        --local-cores 1\
        -pr\
        --cluster "qsub -P {cluster.P}\
                    -l ncpus={cluster.ncpus} \
                    -q {cluster.queue} \
                    -l mem={cluster.mem} \
                    -l wd\
                    -l walltime={cluster.walltime}\
                    -e {cluster.error_out_dir} \
                    -o {cluster.std1_out_dir}" \
        --cluster-config /short/kv78/cellcycle_workflow/cluster.json\
        --rerun-incomplete


