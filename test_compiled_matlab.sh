#!/usr/bin/env bash

export PATH=$(pwd)/matlab/bin:$PATH

run_spm12.sh \
    /home/rogersbp/MATLAB/R2023a \
    function matlab_entrypoint \
    fmriprep1_dir INPUTS/fmriprep-ert1 \
    fmriprep2_dir INPUTS/fmriprep-ert2 \
    hpf_sec 300 \
    psydat_csv OUTPUTS/ert1.csv \
    out_dir OUTPUTS
