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
        --fmriprep_ert1_dir) export fmriprep_ert1_dir="${2}"; shift; shift ;;
        --fmriprep_ert2_dir) export fmriprep_ert2_dir="${2}"; shift; shift ;;
        --fmriprep_rift1_dir) export fmriprep_rift1_dir="${2}"; shift; shift ;;
        --fmriprep_rift2_dir) export fmriprep_rift2_dir="${2}"; shift; shift ;;
        --fmriprep_rift3_dir) export fmriprep_rift3_dir="${2}"; shift; shift ;;
        --fmriprep_rift4_dir) export fmriprep_rift4_dir="${2}"; shift; shift ;;
        --ert_psydat_csv) export ert_psydat_csv="${2}"; shift; shift ;;
        --rift_psydat_csv) export rift_psydat_csv="${2}"; shift; shift ;;
        --hpf_sec) export hpf_sec="${2}"; shift; shift ;;
        --fwhm_mm) export fwhm_mm="${2}"; shift; shift ;;
        --out_dir) export out_dir="${2}"; shift; shift ;;
        *) echo "Input ${1} not recognized" ; shift ;;
    esac
done

# Run the matlab part of the pipeline in xvfb
xvfb-run -n $(($$ + 99)) -s '-screen 0 1600x1200x24 -ac +extension GLX' \
    run_spm12.sh ${MATLAB_RUNTIME} function matlab_entrypoint \
    fmriprep_ert1_dir "${fmriprep_ert1_dir}" \
    fmriprep_ert2_dir "${fmriprep_ert2_dir}" \
    fmriprep_rift1_dir "${fmriprep_rift1_dir}" \
    fmriprep_rift2_dir "${fmriprep_rift2_dir}" \
    fmriprep_rift3_dir "${fmriprep_rift3_dir}" \
    fmriprep_rift4_dir "${fmriprep_rift4_dir}" \
    ert_psydat_csv "${ert_psydat_csv}" \
    rift_psydat_csv "${rift_psydat_csv}" \
    hpf_sec "${hpf_sec}" \
    fwhm_mm "${fwhm_mm}" \
    out_dir "${out_dir}"

# Create QC PDF
convert \
    "${out_dir}"/first_level_design_ert_001.png \
    "${out_dir}"/first_level_result_ert_001.png \
    "${out_dir}"/first_level_report_ert.pdf

# Zip nii
gzip "${out_dir}"/spm_ert/*.nii
