#!/usr/bin/env bash
#
# Primary entrypoint

echo Running $(basename "${BASH_SOURCE}")

# Initialize defaults
export hpf_sec=300
export fwhm_mm=6
export out_dir=/OUTPUTS

# Parse input options
while [[ $# -gt 0 ]]; do
    key="${1}"
    case $key in   
        --fmriprep1_dir) export fmriprep1_dir="${2}"; shift; shift ;;
        --fmriprep2_dir) export fmriprep2_dir="${2}"; shift; shift ;;
        --psydat_file) export psydat_file="${2}"; shift; shift ;;
        --hpf_sec) export hpf_sec="${2}"; shift; shift ;;
        --fwhm_mm) export fwhm_mm="${2}"; shift; shift ;;
        --out_dir) export out_dir="${2}"; shift; shift ;;
        *) echo "Input ${1} not recognized" ; shift ;;
    esac
done

# Convert the psydat file to csv
convert-psydat.py --psydat_file "${psydat_file}" --out_file "${out_dir}"/taskinfo.csv

# Run the matlab part of the pipeline in xvfb
xvfb-run -n $(($$ + 99)) -s '-screen 0 1600x1200x24 -ac +extension GLX' \
    run_spm12.sh function matlab_entrypoint \
    fmriprep1_dir "${fmriprep1_dir}" \
    fmriprep2_dir "${fmriprep2_dir}" \
    psydat_csv "${out_dir}"/taskinfo.csv \
    hpf_sec "${hpf_sec}" \
    fwhm_mm "${fwhm_mm}" \
    out_dir "${out_dir}"

# Create QC PDF
