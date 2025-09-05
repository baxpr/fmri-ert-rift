#!/usr/bin/env bash

docker run \
    --mount type=bind,src=$(pwd -P)/INPUTS,dst=/INPUTS \
    --mount type=bind,src=$(pwd -P)/OUTPUTS,dst=/OUTPUTS \
    fmri-ert-rift:test \
	--fmriprep1_dir /INPUTS/fmriprep-ert1 \
	--fmriprep2_dir /INPUTS/fmriprep-ert2 \
	--hpf_sec 300 \
	--psydat_file /INPUTS/ert1.csv \
	--out_dir /OUTPUTS
