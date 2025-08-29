#!/usr/bin/env bash

export PATH=$(pwd)/bin:$PATH

run_spm12.sh \
    /home/rogersbp/MATLAB/R2023a \
    function matlab_entrypoint \
    fmriprep1_dir $(pwd)/../INPUTS/fmriprep-ert1 \
    fmriprep2_dir $(pwd)/../INPUTS/fmriprep-ert2 \
    hpf_sec 300 \
    psydat_csv $(pwd)/../OUTPUTS/ert1.csv \
    out_dir $(pwd)/../OUTPUTS
