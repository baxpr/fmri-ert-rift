#!/usr/bin/env bash

docker run \
    --mount type=bind,src=$(pwd -P)/INPUTS,dst=/INPUTS \
    --mount type=bind,src=$(pwd -P)/OUTPUTS,dst=/OUTPUTS \
    fmri-ert-rift:test \
    --fmriprep_ert1_dir /INPUTS/fmriprep-ert1/fmriprepBIDS \
    --fmriprep_ert2_dir /INPUTS/fmriprep-ert2/fmriprepBIDS \
    --fmriprep_rift1_dir /INPUTS/fmriprep-rift1/fmriprepBIDS \
    --fmriprep_rift2_dir /INPUTS/fmriprep-rift2/fmriprepBIDS \
    --fmriprep_rift3_dir /INPUTS/fmriprep-rift3/fmriprepBIDS \
    --fmriprep_rift4_dir /INPUTS/fmriprep-rift4/fmriprepBIDS \
    --ert_psydat_csv /OUTPUTS/ert.csv \
    --rift_psydat_csv /OUTPUTS/rift.csv \
	--hpf_sec 300 \
    --fwhm 6 \
	--out_dir /OUTPUTS
