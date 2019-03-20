#!/bin/bash
__conda_setup="$('/home/skurscheid/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/skurscheid/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/skurscheid/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/skurscheid/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

~/miniconda3/envs/snakemake-testing/bin/snakemake -s rules/run_macs2.smk all\
	--configfile ./config.json \
	--use-conda\
	--cluster "qsub -pe threads {cluster.threads} \
                      -q {cluster.queue} \
                      -l virtual_free={cluster.virtual_free} \
                      -l h_vmem={cluster.h_vmem} \
                      -e {cluster.error_out_dir} \
                      -o {cluster.std1_out_dir}" \
	--jobs 32\
	 -d /home/skurscheid/data/experiment\
	--rerun-incomplete\
	--cluster-config /home/skurscheid/workflows/cellcycle_workflow/cluster.json\
	-pr\
