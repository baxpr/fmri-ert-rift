function matlab_entrypoint(varargin)

%% Just quit, if requested - needed for docker build
if numel(varargin)==1 && strcmp(varargin{1},'quit') && isdeployed
	disp('Exiting as requested')
	exit
end
    
% Parse inputs
P = inputParser;
addOptional(P,'fmri1_nii','/OUTPUTS/ctrrfmri1.nii')
addOptional(P,'fmri2_nii','/OUTPUTS/ctrrfmri2.nii')
addOptional(P,'motpar1_txt','/OUTPUTS/rfmri1.par');
addOptional(P,'motpar2_txt','/OUTPUTS/rfmri2.par');
addOptional(P,'meanfmri_nii','/OUTPUTS/ctrrfmri_mean_all.nii')
addOptional(P,'deffwd_nii','/OUTPUTS/y_deffwd.nii')
addOptional(P,'refimg_nii','avg152T1.nii')
addOptional(P,'biasnorm_nii','/OUTPUTS/biasnorm.nii')
addOptional(P,'trialreport_csv','/OUTPUTS/trialreport.csv')
addOptional(P,'hpf_sec','200')
addOptional(P,'fwhm_mm','6')
addOptional(P,'out_dir','/OUTPUTS');
parse(P,varargin{:});
disp(P.Results)

% Run the actual pipeline
matlab_main(P.Results);

% Exit
if isdeployed
	exit
end
