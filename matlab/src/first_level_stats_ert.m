function first_level_stats_ert(inp)

% Input variables
%
%   psydat_csv
%   out_dir
%   fmriprep1_dir
%   fmriprep2_dir
%   hpf_sec

%   Trial phases are Instruction, Image, Response (3)
%   Additionally, model trial type (6)

%    image_category      strategy  
%    ______________    ____________
%
%     {'negative'}     {'ACCEPT'  }
%     {'negative'}     {'AVOID'   }
%     {'negative'}     {'DISTRACT'}
%     {'negative'}     {'LOOK'    }
%     {'negative'}     {'REFRAME' }
%     {'neutral' }     {'LOOK'    }


%% Get timing info from converted psydat
timings = readtable(inp.psydat_csv);

% Add condition column
timings.condition = strcat(timings.image_category,'_',timings.strategy);

% Condition and timing variables
condvars = {'condition'};
timevars = { ...
    'instructional_cue_started', ...
    'instructional_cue_stopped', ...
    'image_started', ...
    'image_stopped', ...
    'affect_rating_ert_started', ...
    'affect_rating_ert_stopped' ...
    };
keepvars = [condvars timevars];

% Find scan start times
scanstarts = sort(timings.startedScanning(~isnan(timings.startedScanning)));

% Run-specific timing info relative to beginning of first fmri
clear trialtimes
for r = 1:2
    blockfile = ['ert_block_' num2str(r) '.csv'];
    trialtimes{r} = timings(strcmp(timings.block_file,blockfile),keepvars);
    trialtimes{r}(:,timevars) = trialtimes{r}(:,timevars) - scanstarts(r);
end


%% fmriprep stuff

% Scale motion params and save in SPM friendly format
for r = 1:2
    confD = dir([inp.(['fmriprep' num2str(r) '_dir']) '/sub*/ses*/func/*_desc-confounds_timeseries.tsv']);
    conf = readtable(fullfile(confD(1).folder,confD(1).name),'FileType','text','Delimiter','tab');
    motT = conf(:,{'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'});
    mot = zscore(table2array(motT));
    writematrix(mot, fullfile(inp.out_dir,['motpar' num2str(r) '.txt']))
end

% Find preprocessed image files
clear niigz
for r = 1:2
    niigzD = dir([inp.(['fmriprep' num2str(r) '_dir']) '/sub*/ses*/func/*_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz']);
    niigz{r} = fullfile(niigzD(1).folder,niigzD(1).name);
end


%% Other params needed by SPM

% Filter param
hpf_sec = str2double(inp.hpf_sec);

% Get TRs and check
N = nifti(inp.swfmri1_nii);
tr = N.timing.tspace;
for r = 2
	N = nifti(inp.(['swfmri' num2str(r) '_nii']));
	if abs(N.timing.tspace-tr) > 0.001
		error('TR not matching for run %d',r)
	end
end



%% OLD BELOW HERE



%% Design
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = ...
	{fullfile(inp.out_dir,['spm_' tag])};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
matlabbatch{1}.spm.stats.fmri_spec.mask = {[spm('dir') '/tpm/mask_ICV.nii']};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

for r = 1:2
	
	% Initialize
	thist = trials(trials.Run==r,:);
	ind = ismember(thist.Outcome,{'Win','Lose'});
	c = 0;
	
	% Session-specific scans, regressors, params
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = ...
		{inp.(['swfmri' num2str(r) '_nii'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = ...
		struct('name', {}, 'val', {});
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = ...
		{fullfile(inp.out_dir,['motpar' num2str(r) '.txt'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf_sec;
	
	% Condition: Cue
	c = c + 1;
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).name = 'Cue';
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).onset = ...
		thist.T1_TrialStart_fMRIsec(ind);
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).duration = 0;
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).tmod = 0;

	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod(1).name = 'mu33';
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod(1).param = thist.traj_mu_33(ind);
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod(1).poly = 1;

	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).orth = 1;

	% Condition: Feedback
	c = c + 1;
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).name = 'Feedback';
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).onset = ...
		thist.T3_FeedbackOnset_fMRIsec(ind);
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).duration = 0;
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).tmod = 0;

	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod(1).name = 'epsi3';
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod(1).param = thist.traj_epsi_3(ind);
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod(1).poly = 1;

	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).orth = 1;	
	
end

%% Estimate
matlabbatch{2}.spm.stats.fmri_est.spmmat = ...
	fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


%% Contrasts
%
% Parameters per session SPM.xX.name' are
%    {'Sn(1) Cue*bf(1)'             }
%    {'Sn(1) Cuexmu33^1*bf(1)'      }
%    {'Sn(1) Feedback*bf(1)'        }
%    {'Sn(1) Feedbackxepsi3^1*bf(1)'}
%    {'Sn(1) R1'                    }
%    {'Sn(1) R2'                    }
%    {'Sn(1) R3'                    }
%    {'Sn(1) R4'                    }
%    {'Sn(1) R5'                    }
%    {'Sn(1) R6'                    }
matlabbatch{3}.spm.stats.con.spmmat = ...
	matlabbatch{2}.spm.stats.fmri_est.spmmat;
matlabbatch{3}.spm.stats.con.delete = 1;
c = 0;

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Cue';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Cue x Mu33';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 1 0 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Feedback';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Feedback x Epsi3';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Cue + Feedback';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0.5 0 0.5 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';


% Inverse of all existing contrasts since SPM won't show us both sides
numc = numel(matlabbatch{3}.spm.stats.con.consess);
for k = 1:numc
        c = c + 1;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = ...
                ['Neg ' matlabbatch{3}.spm.stats.con.consess{c-numc}.tcon.name];
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
                - matlabbatch{3}.spm.stats.con.consess{c-numc}.tcon.weights;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';
end


%% Review design
matlabbatch{4}.spm.stats.review.spmmat = ...
	matlabbatch{2}.spm.stats.fmri_est.spmmat;
matlabbatch{4}.spm.stats.review.display.matrix = 1;
matlabbatch{4}.spm.stats.review.print = false;

matlabbatch{5}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = ...
        fullfile(inp.out_dir,['first_level_design_' tag '.png']);
matlabbatch{5}.cfg_basicio.run_ops.call_matlab.outputs = cell(1,0);
matlabbatch{5}.cfg_basicio.run_ops.call_matlab.fun = 'spm_window_print';


%% Save batch and run
save(fullfile(inp.out_dir,['spmbatch_first_level_stats_' tag '.mat']),'matlabbatch')
spm_jobman('run',matlabbatch);

% And save contrast names
numc = numel(matlabbatch{3}.spm.stats.con.consess);
connames = table((1:numc)','VariableNames',{'ConNum'});
for k = 1:numc
	try
		connames.ConName{k,1} = ...
			matlabbatch{3}.spm.stats.con.consess{k}.tcon.name;
	catch
		connames.ConName{k,1} = ...
			matlabbatch{3}.spm.stats.con.consess{k}.fcon.name;
	end
end
writetable(connames,fullfile(inp.out_dir,['spm_contrast_names_' tag '.csv']));


%% Results display
% Needed to create the spmT even if we don't get the figure window
xSPM = struct( ...
    'swd', matlabbatch{1}.spm.stats.fmri_spec.dir, ...
    'title', '', ...
    'Ic', 1, ...
    'n', 0, ...
    'Im', [], ...
    'pm', [], ...
    'Ex', [], ...
    'u', 0.005, ...
    'k', 10, ...
    'thresDesc', 'none' ...
    );
[hReg,xSPM] = spm_results_ui('Setup',xSPM);

% Show on the subject MNI anat
spm_sections(xSPM,hReg,inp.biasnorm_nii)

% Jump to global max activation
spm_mip_ui('Jump',spm_mip_ui('FindMIPax'),'glmax');

