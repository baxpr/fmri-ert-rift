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

% SPM init
spm_jobman('initcfg');

% Run the actual pipeline
outp = fileprep_ert_rift(P.Results);
first_level_stats_ert(outp);

% Exit
if isdeployed
	exit
end
