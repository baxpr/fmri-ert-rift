#!/usr/bin/env bash

export PATH=$(pwd)/bin:$PATH

run_matlab_entrypoint.sh \
    /home/rogersbp/MATLAB/R2023a \
	fmriprep1_dir $(pwd)/INPUTS/fmriprep-ert1 \
	fmriprep2_dir $(pwd)/INPUTS/fmriprep-ert2 \
	hpf_sec 300 \
	psydat_csv $(pwd)/OUTPUTS/ert1.csv \
	out_dir $(pwd)/OUTPUTS
