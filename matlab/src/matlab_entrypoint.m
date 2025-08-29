function matlab_entrypoint(varargin)

%% Just quit, if requested - needed for docker build
if numel(varargin)==1 && strcmp(varargin{1},'quit') && isdeployed
	disp('Exiting as requested')
	exit
end
    
% Parse inputs
P = inputParser;
addOptional(P,'fmriprep1_dir','/OUTPUTS/fmriprep1')
addOptional(P,'fmriprep2_dir','/OUTPUTS/fmriprep2')
addOptional(P,'psydat_csv','/OUTPUTS/psydat.csv')
addOptional(P,'hpf_sec','300')
addOptional(P,'fwhm_mm','6')
addOptional(P,'out_dir','/OUTPUTS');
parse(P,varargin{:});
disp(P.Results)

% SPM init
spm_jobman('initcfg');

% Run the actual pipeline
first_level_stats_ert(P.Results);

% Exit
if isdeployed
	exit
end
