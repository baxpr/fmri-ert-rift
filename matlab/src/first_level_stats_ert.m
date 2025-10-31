function first_level_stats_ert(inp)


%% Other params needed by SPM

spm('defaults','fmri')

% Filter param
hpf_sec = str2double(inp.hpf_sec);

% Get TRs and check across runs
N = nifti(inp.fmri_nii{1});
tr = N.timing.tspace;
for r = 2:numel(inp.fmri_nii)
    N = nifti(inp.fmri_nii{r});
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

% Smooth fmriprep's fmri timeseries and get smoothed filenames
fwhm_mm = str2double(inp.fwhm_mm);
clear smfri_nii
c = 0;
for imgs = inp.fmri_nii

    clear matlabbatch
    matlabbatch{1}.spm.spatial.smooth.data = imgs(1);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm_mm fwhm_mm fwhm_mm];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',matlabbatch);

    [~,n,e] = fileparts(imgs{1});
    c = c + 1;
    sfmri_nii{c} = fullfile(inp.out_dir,['s' n e]);

end


% General design
clear matlabbatch
spmdir = fullfile(inp.out_dir,'spm_ert');
matlabbatch{1}.spm.stats.fmri_spec.dir = {spmdir};
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


for r = 1:6

    % Session-specific scans, regressors, params
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = sfmri_nii(r);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = inp.motpar_txt(r);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf_sec;

    % Conditions per run
    k = 0;
    for c = 1:numel(conds)

        % Only keep specific condition, ERT-type trials
        inds = strcmp(inp.trialtimes{r}.ACTUAL_CONDITION,conds{c}) ...
            & strcmp(inp.trialtimes{r}.TRIAL_TYPE_ERT,'ERT');
        %fprintf('Found %d %s trials run %d\n',sum(inds),conds{c},r);
        if sum(inds)==0, continue, end

        % Condition image
        k = k + 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).name = [conds{c} '_Image'];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).onset = ...
            inp.trialtimes{r}.image_started(inds);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).duration = ...
            inp.trialtimes{r}.image_stopped(inds) ...
            - inp.trialtimes{r}.image_started(inds);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).orth = 1;

        % Condition response
        k = k + 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).name = [conds{c} '_Response'];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).onset = ...
            inp.trialtimes{r}.affect_rating_started(inds);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).duration = ...
            inp.trialtimes{r}.affect_rating_stopped(inds) ...
            - inp.trialtimes{r}.affect_rating_started(inds);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).pmod = ...
            struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(k).orth = 1;

    end

end


%% Estimate
spmmat = fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.spmmat = spmmat;
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


%% Run, to get info we need for the contrasts
save(fullfile(inp.out_dir,'spmbatch_first_level_stats_ert_spec_est.mat'),'matlabbatch')
spm_jobman('run',matlabbatch);


%% Contrasts
%
% Parameters per session SPM.xX.name vary.
%
% Desired contrasts are "look_NEGATIVE" vs each of the other five

% Get predictor names (SPM.xX.name contains 'negative_LOOK_Image' etc)
load(spmmat{1},'SPM');

% We will weight contrast entries by number of trials in the session

% Find trial counts using SPM.xX.name
trialcounts = table();
t = 0;
for p = 1:numel(SPM.xX.name)
    t = t + 1;
    trialcounts.predictornum(t,1) = p;
    trialcounts.predictorname(t,1) = SPM.xX.name(p);
    sn = regexp(SPM.xX.name{p},'Sn\((.)\)','tokens');
    trialcounts.session(t,1) = str2num(sn{1}{1});
    cn = regexp(SPM.xX.name{p},'Sn\(.\) (.*)\*','tokens');
    if numel(cn)>0
        trialcounts.condition{t,1} = cn{1}{1};
    else
        trialcounts.condition{t,1} = '';
    end
    u = find(strcmp([SPM.Sess(trialcounts.session(t,1)).U.name],trialcounts.condition{t,1}));
    if numel(u)>0
        trialcounts.trialcount(t,1) = numel(SPM.Sess(trialcounts.session(t,1)).U(u).ons);
    else
        trialcounts.trialcount(t,1) = 0;
    end
end

% Setup
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = spmmat;
matlabbatch{1}.spm.stats.con.delete = 1;
c = 0;

% Sanity check, Image vs Response
c = c + 1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = 'Image gt Response';
matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = zeros(1,numel(trialcounts.predictorname));
ind1 = find(contains(trialcounts.condition,'_Image'));
cnt1 = trialcounts.trialcount(ind1);
wt1 = cnt1 ./ sum(cnt1);
ind2 = find(contains(trialcounts.condition,'_Response'));
cnt2 = trialcounts.trialcount(ind2);
wt2 = cnt2 ./ sum(cnt2);
matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind1) = wt1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind2) = -wt2;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';

% Individual conditions
conds = {
    'negative_LOOK'
    'neutral_LOOK'
    'negative_ACCEPT'
    'negative_AVOID'
    'negative_DISTRACT'
    'negative_REFRAME'
    };
types = {'Image', 'Response'};
for a = 1:numel(conds)
    for b = 1:numel(types)
        c = c + 1;
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = [conds{a} '_' types{b}];
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = zeros(1,numel(trialcounts.predictorname));
        ind1 = find(contains(trialcounts.condition,[conds{a} '_' types{b}]));
        cnt1 = trialcounts.trialcount(ind1);
        wt1 = cnt1 ./ sum(cnt1);
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind1) = wt1;
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
    end
end

% Each condition vs negative_LOOK for both Image and Response
conds = {
    'neutral_LOOK'
    'negative_ACCEPT'
    'negative_AVOID'
    'negative_DISTRACT'
    'negative_REFRAME'
    };
types = {'Image', 'Response'};
for a = 1:numel(conds)
    for b = 1:numel(types)
        c = c + 1;
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = [conds{a} '_gt_negative_LOOK_' types{b}];
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = zeros(1,numel(trialcounts.predictorname));
        ind1 = find(contains(trialcounts.condition,[conds{a} '_' types{b}]));
        cnt1 = trialcounts.trialcount(ind1);
        wt1 = cnt1 ./ sum(cnt1);
        ind2 = find(contains(trialcounts.condition,['negative_LOOK_' types{b}]));
        cnt2 = trialcounts.trialcount(ind2);
        wt2 = cnt2 ./ sum(cnt2);
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind1) = wt1;
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind2) = -wt2;
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
    end
end

% accept + reframe (vs negative) and distract + avoidance (vs negative) 
types = {'Image', 'Response'};
for b = 1:numel(types)
    c = c + 1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = ['AVOID_DISTRACT_gt_negative_LOOK_' types{b}];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = zeros(1,numel(trialcounts.predictorname));
    ind1 = find(contains(trialcounts.condition,['negative_AVOID_' types{b}]));
    cnt1 = trialcounts.trialcount(ind1);
    wt1 = cnt1 ./ sum(cnt1);
    ind2 = find(contains(trialcounts.condition,['negative_DISTRACT_' types{b}]));
    cnt2 = trialcounts.trialcount(ind2);
    wt2 = cnt2 ./ sum(cnt2);
    ind3 = find(contains(trialcounts.condition,['negative_LOOK_' types{b}]));
    cnt3 = trialcounts.trialcount(ind3);
    wt3 = cnt3 ./ sum(cnt3);
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind1) = 0.5 * wt1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind2) = 0.5 * wt2;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind3) = -wt3;
end

types = {'Image', 'Response'};
for b = 1:numel(types)
    c = c + 1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = ['ACCEPT_REFRAME_gt_negative_LOOK_' types{b}];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = zeros(1,numel(trialcounts.predictorname));
    ind1 = find(contains(trialcounts.condition,['negative_ACCEPT_' types{b}]));
    cnt1 = trialcounts.trialcount(ind1);
    wt1 = cnt1 ./ sum(cnt1);
    ind2 = find(contains(trialcounts.condition,['negative_REFRAME_' types{b}]));
    cnt2 = trialcounts.trialcount(ind2);
    wt2 = cnt2 ./ sum(cnt2);
    ind3 = find(contains(trialcounts.condition,['negative_LOOK_' types{b}]));
    cnt3 = trialcounts.trialcount(ind3);
    wt3 = cnt3 ./ sum(cnt3);
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind1) = 0.5 * wt1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind2) = 0.5 * wt2;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind3) = -wt3;
end

% accept + reframe  vs  distract + avoidance
types = {'Image', 'Response'};
for b = 1:numel(types)
    c = c + 1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = ['ACCEPT_REFRAME_gt_AVOID_DISTRACT_' types{b}];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = zeros(1,numel(trialcounts.predictorname));
    ind1 = find(contains(trialcounts.condition,['negative_ACCEPT_' types{b}]));
    cnt1 = trialcounts.trialcount(ind1);
    wt1 = cnt1 ./ sum(cnt1);
    ind2 = find(contains(trialcounts.condition,['negative_REFRAME_' types{b}]));
    cnt2 = trialcounts.trialcount(ind2);
    wt2 = cnt2 ./ sum(cnt2);
    ind3 = find(contains(trialcounts.condition,['negative_AVOID_' types{b}]));
    cnt3 = trialcounts.trialcount(ind3);
    wt3 = cnt3 ./ sum(cnt3);
    ind4 = find(contains(trialcounts.condition,['negative_DISTRACT_' types{b}]));
    cnt4 = trialcounts.trialcount(ind4);
    wt4 = cnt4 ./ sum(cnt4);
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind1) = 0.5 * wt1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind2) = 0.5 * wt2;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind3) = -0.5 * wt3;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights(ind4) = -0.5 * wt4;
end


%% Inverse of all existing contrasts since SPM won't show us both sides
numc = numel(matlabbatch{1}.spm.stats.con.consess);
for k = 1:numc
    c = c + 1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = ...
        ['Neg ' matlabbatch{1}.spm.stats.con.consess{c-numc}.tcon.name];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = ...
        - matlabbatch{1}.spm.stats.con.consess{c-numc}.tcon.weights;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
end


%% Review design
matlabbatch{2}.spm.stats.review.spmmat = spmmat;
matlabbatch{2}.spm.stats.review.display.matrix = 1;
matlabbatch{2}.spm.stats.review.print = false;

matlabbatch{3}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = ...
    fullfile(inp.out_dir,'first_level_design_ert.png');
matlabbatch{3}.cfg_basicio.run_ops.call_matlab.outputs = cell(1,0);
matlabbatch{3}.cfg_basicio.run_ops.call_matlab.fun = 'spm_window_print';


%% Save batch and run
save(fullfile(inp.out_dir,'spmbatch_first_level_stats_ert_con_run_review.mat'),'matlabbatch')
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
    'swd', spmdir, ...
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
spm_sections(xSPM,hReg,inp.atlasT1_nii)

% Jump to activation location
%spm_mip_ui('Jump',spm_mip_ui('FindMIPax'),'glmax');  % Global max
spm_mip_ui('SetCoords',[5 -96 6]);  % Task specific MNI location

% Screenshot
spm_window_print(fullfile(inp.out_dir,'first_level_result_ert.png'));

