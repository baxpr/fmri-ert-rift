function matlab_entrypoint(varargin)

%% Just quit, if requested - needed for docker build
if numel(varargin)==1 && strcmp(varargin{1},'quit') && isdeployed
	disp('Exiting as requested')
	exit
end
    
% Parse inputs
P = inputParser;
addOptional(P,'fmriprep_ert1_dir','none')
addOptional(P,'fmriprep_ert2_dir','none')
addOptional(P,'fmriprep_rift1_dir','none')
addOptional(P,'fmriprep_rift2_dir','none')
addOptional(P,'fmriprep_rift3_dir','none')
addOptional(P,'fmriprep_rift4_dir','none')
addOptional(P,'ert_psydat_csv','/OUTPUTS/ert1.csv')
addOptional(P,'rift_psydat_csv','/OUTPUTS/rift1.csv')
addOptional(P,'hpf_sec','300')
addOptional(P,'fwhm_mm','6')
addOptional(P,'out_dir','/OUTPUTS');
parse(P,varargin{:});
disp(P.Results)

% FIXME Due to small number of trials and inconsistent presence across
% runs, let's concatenate and treat the full data set as a single "run".
%
% 1. Drop short non-ERT/RIFT trial segments at end of ERT2 and RIFT4 (?)
%
% 2. Account for missing runs ('none' supplied as argument) 
%
% 3. Adjust event timings by adding the number of actual TRs previous in
% the combined run at the beginning of the segment.
%
% 4. Combine fmriprep masks (? or just use SPM TPM as with analysis) and
% scale fmri grand mean per run
%
% 5. Fix the filtering - can't be done overall, needs to be per run then
% turned off for fitting. Hm. Even the run mean can vary at a single voxel,
% even after global grand mean scaling - would need to scale per voxel?
%
%
% Alternative is to build contrasts manually, leaving out any runs that
% don't have events. And possibly leaving out regressors for a single run
% when events don't exist.


% SPM init
spm_jobman('initcfg');

% Run the actual pipeline
first_level_stats_ert(P.Results);

% Exit
if isdeployed
	exit
end
