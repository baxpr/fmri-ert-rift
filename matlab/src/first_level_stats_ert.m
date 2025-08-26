function first_level_stats_ert(inp)

% Input variables
%   psydat_csv
%   fmriprep1_dir
%   fmriprep2_dir
%   hpf_sec
%   out_dir

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
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
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
nruns = 2;
for r = 1:nruns
    blockfile = ['ert_block_' num2str(r) '.csv'];
    trialtimes{r} = timings(strcmp(timings.block_file,blockfile),keepvars);
    trialtimes{r}(:,timevars) = trialtimes{r}(:,timevars) - scanstarts(r);
end


%% Find fmriprep files

% Scale motion params and save in SPM friendly format
for r = 1:nruns
    confD = dir([inp.(['fmriprep' num2str(r) '_dir']) '/sub*/ses*/func/*_desc-confounds_timeseries.tsv']);
    conf = readtable(fullfile(confD(1).folder,confD(1).name),'FileType','text','Delimiter','tab');
    motT = conf(:,{'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'});
    mot = zscore(table2array(motT));
    writematrix(mot, fullfile(inp.out_dir,['motpar' num2str(r) '.txt']))
end

% Find preprocessed image files, copy, and unzip
clear fmri_nii
for r = 1:nruns
    niigzD = dir([inp.(['fmriprep' num2str(r) '_dir']) '/sub*/ses*/func/*_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz']);
    fmri_nii{r} = fullfile(inp.out_dir,['fmri' num2str(r) '.nii']);
    copyfile( ...
        fullfile(niigzD(1).folder,niigzD(1).name), ...
        [fmri_nii{r} '.gz'] ...
        );
end
gunzip(fullfile(inp.out_dir,'fmri*.nii.gz'));
delete(fullfile(inp.out_dir,'fmri*.nii.gz'));

% Also the T1
niigzD = dir([inp.('fmriprep1_dir') '/sub*/ses*/anat/*_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz']);
atlasT1_nii = fullfile(inp.out_dir,'t1.nii');
copyfile( ...
    fullfile(niigzD(1).folder,niigzD(1).name), ...
    [atlasT1_nii '.gz'] ...
    );
gunzip(fullfile(inp.out_dir,'t1.nii.gz'));
delete(fullfile(inp.out_dir,'t1.nii.gz'));

  

%% Other params needed by SPM

% Filter param
hpf_sec = str2double(inp.hpf_sec);

% Get TRs and check across runs
N = nifti(fmri_nii{1});
tr = N.timing.tspace;
for r = 2:nruns
	N = nifti(fmri_nii{r});
	if abs(N.timing.tspace-tr) > 0.001
		error('TR not matching for run %d',r)
	end
end


%% Design

% Conditions pre-specified
conds = {
    'negative_LOOK'
    'neutral_LOOK'
    'negative_ACCEPT'
    'negative_AVOID'
    'negative_DISTRACT'
    'negative_REFRAME'
    };

% General design
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = ...
	{fullfile(inp.out_dir,'spm_ert')};
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


for r = 1:nruns

	% Session-specific scans, regressors, params
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = ...
		fmri_nii(r);
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = ...
		struct('name', {}, 'val', {});
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = ...
		{fullfile(inp.out_dir,['motpar' num2str(r) '.txt'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf_sec;
	
    % Conditions per run
    k = 0;
    for c = 1:numel(conds)

        inds = strcmp(trialtimes{r}.condition,conds{c});

%        % Condition cue
%        k = k + 1;
%    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).name = [conds{c} '_Cue'];
%    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).onset = ...
%    		trialtimes{r}.instructional_cue_started(inds);
%    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).duration = ...
%            trialtimes{r}.instructional_cue_stopped(inds) ...
%            - trialtimes{r}.instructional_cue_started(inds);
%    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).tmod = 0;
%        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).pmod = ...
%            struct('name', {}, 'param', {}, 'poly', {});
%    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).orth = 1;

        % Condition image
        k = k + 1;
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).name = [conds{c} '_Image'];
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).onset = ...
    		trialtimes{r}.image_started(inds);
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).duration = ...
            trialtimes{r}.image_stopped(inds) ...
            - trialtimes{r}.image_started(inds);
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).orth = 1;

        % Condition response
        k = k + 1;
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).name = [conds{c} '_Response'];
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).onset = ...
    		trialtimes{r}.affect_rating_ert_started(inds);
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).duration = ...
            trialtimes{r}.affect_rating_ert_stopped(inds) ...
            - trialtimes{r}.affect_rating_ert_started(inds);
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).orth = 1;

    end

end


%% Estimate
matlabbatch{2}.spm.stats.fmri_est.spmmat = ...
	fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


%% Contrasts
%
% Parameters per session SPM.xX.name' are
%    {'negative_LOOK_Image'       }
%    {'negative_LOOK_Response'    }
%    {'neutral_LOOK_Image'        }
%    {'neutral_LOOK_Response'     }
%    {'negative_ACCEPT_Image'     }
%    {'negative_ACCEPT_Response'  }
%    {'negative_AVOID_Image'      }
%    {'negative_AVOID_Response'   }
%    {'negative_DISTRACT_Image'   }
%    {'negative_DISTRACT_Response'}
%    {'negative_REFRAME_Image'    }
%    {'negative_REFRAME_Response' }
%    {'Sn(1) R1'                  }
%    {'Sn(1) R2'                  }
%    {'Sn(1) R3'                  }
%    {'Sn(1) R4'                  }
%    {'Sn(1) R5'                  }
%    {'Sn(1) R6'                  }
matlabbatch{3}.spm.stats.con.spmmat = ...
	matlabbatch{2}.spm.stats.fmri_est.spmmat;
matlabbatch{3}.spm.stats.con.delete = 1;
c = 0;

% Sanity check
c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Image gt Response';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    1/6 * [1 -1   1 -1   1 -1   1 -1   1 -1   1 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

% Individual conditions
types = {'Image', 'Response'};
for a = 1:numel(conds)
    for b = 1:numel(types)
        c = c + 1;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = [conds{a} '_' types{b}];
        v1 = [zeros(1,b-1) 1 zeros(1,numel(types)-b)];
        cvec = zeros(1,numel(conds)*numel(types));
        ind = a*numel(types)-numel(types)+1;
        cvec(ind:ind+numel(types)-1) = v1;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = cvec;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';
    end
end

% negative_LOOK vs others
c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Image negative_LOOK gt neutral_LOOK';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [1 0   -1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Response negative_LOOK gt neutral_LOOK';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [0 1   0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Image negative_LOOK gt negative_ACCEPT';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [1 0   0 0   -1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Response negative_LOOK gt negative_ACCEPT';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [0 1   0 0   0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Image negative_LOOK gt negative_AVOID';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [1 0   0 0   0 0   -1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Response negative_LOOK gt negative_AVOID';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [0 1   0 0   0 0   0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Image negative_LOOK gt negative_DISTRACT';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [1 0   0 0   0 0   0 0   -1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Response negative_LOOK gt negative_DISTRACT';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [0 1   0 0   0 0   0 0   0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Image negative_LOOK gt negative_REFRAME';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [1 0   0 0   0 0   0 0   0 0   -1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Response negative_LOOK gt negative_REFRAME';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
    [0 1   0 0   0 0   0 0   0 0   0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';


%
%% Inverse of all existing contrasts since SPM won't show us both sides
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
        fullfile(inp.out_dir,'first_level_design_ert.png');
matlabbatch{5}.cfg_basicio.run_ops.call_matlab.outputs = cell(1,0);
matlabbatch{5}.cfg_basicio.run_ops.call_matlab.fun = 'spm_window_print';


%% Save batch and run
save(fullfile(inp.out_dir,'spmbatch_first_level_stats_ert.mat'),'matlabbatch')
spm('defaults','fmri')
spm_jobman('run',matlabbatch);

% And save contrast names
%numc = numel(matlabbatch{3}.spm.stats.con.consess);
%connames = table((1:numc)','VariableNames',{'ConNum'});
%for k = 1:numc
%	try
%		connames.ConName{k,1} = ...
%			matlabbatch{3}.spm.stats.con.consess{k}.tcon.name;
%	catch
%		connames.ConName{k,1} = ...
%			matlabbatch{3}.spm.stats.con.consess{k}.fcon.name;
%	end
%end
%writetable(connames,fullfile(inp.out_dir,['spm_contrast_names_' tag '.csv']));


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
spm_sections(xSPM,hReg,atlasT1_nii)

% Jump to activation location
%spm_mip_ui('Jump',spm_mip_ui('FindMIPax'),'glmax');  % Global max
spm_mip_ui('SetCoords',[5 -96 6]);  % Task specific MNI location

