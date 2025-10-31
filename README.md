# ERT task

## Inputs

- fMRI in atlas space, preprocessed with fmriprep, for the six task runs (ERT_1, ERT_2, RIFT_1, RIFT_2, RIFT_3, RIFT_4)

- CSV format task timing files for the ERT and RIFT sections of the task. These are the output of https://github.com/VUIIS/psychopy-iis


## Methods

For the ERT task, these conditions are modeled:

    Valence      Task
    
    neutral      LOOK
    negative     LOOK
    negative     ACCEPT
    negative     AVOID
    negative     DISTRACT
    negative     REFRAME

For each condition, the image presentation period and the response period are modeled separately.

Each period is modeled as a block convolved with the standard SPM12 hemodynamic response function.

The second trial of every RIFT pair in the four RIFT runs is ignored and not modeled.

The number of trials per run varies, so contrasts are weighted across runs such that each individual trial receives equal weight.

The six estimated motion parameters for each run are included in the model as covariates of no interest.


## Outputs

Several `spmbatch` files containing the design and analysis configuration.

Standard SPM12 results directory, `spm_ert`

A PDF suitable for quick visual QC review, `first_level_report_ert.pdf`.


