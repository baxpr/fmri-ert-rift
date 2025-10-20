function outp = fileprep_ert_rift(inp)

% Input variables
%   ert_psydat_csv
%   rift_psydat_csv
%   fmriprep_ert1_dir
%   fmriprep_ert2_dir
%   fmriprep_rift1_dir
%   fmriprep_rift2_dir
%   fmriprep_rift3_dir
%   fmriprep_rift4_dir
%   hpf_sec
%   fwhm_mm
%   out_dir

% Output variables listed at end


%% ERT - Get timing info from converted psydat
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
timings_ert = readtable(inp.ert_psydat_csv);

% Add condition column
timings_ert.ACTUAL_CONDITION = strcat(timings_ert.image_category,'_',timings_ert.strategy);

% Condition and timing variables
condvars = { ...
    'block_file', ...
    'image_category', ...
    'strategy', ...
    'image_filename', ...
    'ACTUAL_CONDITION' ...
    };
timevars = { ...
    'instructional_cue_started', ...
    'instructional_cue_stopped', ...
    'image_started', ...
    'image_stopped', ...
    'affect_rating_ert_started', ...
    'affect_rating_ert_stopped' ...
    };
keepvars = [condvars timevars];

% Find scan start times. startedScanning has inconsistent/unknown value in
% some cases so don't use it
scanstarts = sort(timings_ert.wait_for_scanner_stopped(~isnan(timings_ert.wait_for_scanner_stopped)));

% Run-specific timing info relative to beginning of first fmri
clear trialtimes_ert
nruns = 2;
for r = 1:nruns
    blockfile = ['ert_block_' num2str(r) '.csv'];
    trialtimes_ert{r} = timings_ert(strcmp(timings_ert.block_file,blockfile),keepvars);
    trialtimes_ert{r}(:,timevars) = trialtimes_ert{r}(:,timevars) - scanstarts(r);
    trialtimes_ert{r}.TRIAL_TYPE_ERT(:) = {'ERT'};

    % Few empties to match the rift version
    trialtimes_ert{r}.new_strategy_chosen(:) = {''};
    trialtimes_ert{r}.switched(:) = {''};
    trialtimes_ert{r}.TRIAL_TYPE_RIFT(:) = {''};
    trialtimes_ert{r}.TRIAL_TYPE_RIFT_PAIR(:) = {''};

    % Tweak names
    trialtimes_ert{r}.affect_rating_started = trialtimes_ert{r}.affect_rating_ert_started;
    trialtimes_ert{r}.affect_rating_stopped = trialtimes_ert{r}.affect_rating_ert_stopped;

    trialtimes_ert{r} = trialtimes_ert{r}(:, { ...
        'block_file', ...
        'image_category', ...
        'strategy', ...
        'new_strategy_chosen', ...
        'image_filename', ...
        'switched', ...
        'TRIAL_TYPE_ERT', ...
        'TRIAL_TYPE_RIFT', ...
        'TRIAL_TYPE_RIFT_PAIR', ...
        'ACTUAL_CONDITION', ...
        'instructional_cue_started', ...
        'instructional_cue_stopped', ...
        'image_started', ...
        'image_stopped', ...
        'affect_rating_started', ...
        'affect_rating_stopped' ...
        });

    writetable(trialtimes_ert{r},fullfile(inp.out_dir,['trialtimes_ERT_run' num2str(r)]));

end


%% RIFT - Get timing info from converted psydat
timings_rift = readtable(inp.rift_psydat_csv);

% Add condition column
timings_rift.LABELED_CONDITION = strcat(timings_rift.image_category,'_',timings_rift.strategy);

% Condition and timing variables
condvars = { ...
    'block_file', ...
    'image_category', ...
    'strategy', ...
    'LABELED_CONDITION', ...
    'switched', ...
    'new_strategy_chosen', ...
    'image_filename' ...
    };
timevars = { ...
    'instructional_cue_rift_started', ...
    'instructional_cue_rift_stopped', ...
    'image_2_started', ...
    'image_2_stopped', ...
    'affect_rating_rift_started', ...
    'affect_rating_rift_stopped' ...
    };
keepvars = [condvars timevars];

% Find scan start times. startedScanning has inconsistent/unknown value in
% some cases so don't use it
scanstarts = sort(timings_rift.wait_for_scanner_stopped(~isnan(timings_rift.wait_for_scanner_stopped)));

% Run-specific timing info relative to beginning of first fmri
clear trialtimes_rift
nruns = 4;
for r = 1:nruns
    blockfile = ['rift_block_' num2str(r) '.csv'];
    trialtimes_rift{r} = timings_rift(strcmp(timings_rift.block_file,blockfile),keepvars);
    trialtimes_rift{r}(:,timevars) = trialtimes_rift{r}(:,timevars) - scanstarts(r);

    % new_strategy_chosen only makes sense if they switched (otherwise is
    % historical)
    % FIXME does the new strat apply to the NEXT trial? I.e. edit those
    % conditions? Or what? Can't tell which trial times are what condition
    trialtimes_rift{r}.new_strategy_chosen(~strcmp(trialtimes_rift{r}.switched,'yes')) = {''};

    % Some trials get two rows, but only one row has timing info
    trialtimes_rift{r} = trialtimes_rift{r}(~isnan(trialtimes_rift{r}.instructional_cue_rift_started),:);

    % if image is identical between N and N+1,
    %    mark N as ERT / "normal" trial
    %    mark N+1 as RIFT trial
    %    if N's switch is 'yes',
    %       also rename N+1 condition by new_strategy_chosen
    %
    %  ???

    trialtimes_rift{r}.shifted_filename(:) = {''};
    trialtimes_rift{r}.shifted_filename(2:end) = trialtimes_rift{r}.image_filename(1:end-1);

    % Default to ERT label
    trialtimes_rift{r}.TRIAL_TYPE_RIFT(:) = {''};

    % Second trial of RIFT pair labeled RIFT2
    trialtimes_rift{r}.TRIAL_TYPE_RIFT(strcmp(trialtimes_rift{r}.image_filename,trialtimes_rift{r}.shifted_filename)) = {'RIFT2'};

    % Anything before a RIFT2 is a RIFT1
    k = find(strcmp(trialtimes_rift{r}.TRIAL_TYPE_RIFT,'RIFT2'));
    trialtimes_rift{r}.TRIAL_TYPE_RIFT(k-1) = {'RIFT1'};

    % Unlabeled and RIFT1 are also used for ERT
    trialtimes_rift{r}.TRIAL_TYPE_ERT(:) = {'ERT'};
    trialtimes_rift{r}.TRIAL_TYPE_ERT(strcmp(trialtimes_rift{r}.TRIAL_TYPE_RIFT,'RIFT2')) = {''};

    % Fix condition labels
    trialtimes_rift{r}.ACTUAL_CONDITION = trialtimes_rift{r}.LABELED_CONDITION;
    update_conds = find(strcmp(trialtimes_rift{r}.switched,'yes')) + 1;
    trialtimes_rift{r}.ACTUAL_CONDITION(update_conds) = ...
        strcat( ...
        trialtimes_rift{r}.image_category(update_conds), ...
        '_', ...
        trialtimes_rift{r}.new_strategy_chosen(update_conds-1) ...
        );

    % Also label RIFT trials as these (which should be all the
    % possibilities):
    %    RtoA_1, RtoA_2
    %    AtoR_1, AtoR_2
    %    RtoR_1, RtoR_2
    %    AtoA_1, AtoA_2
    trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR(:) = {''};
    for k = find(strcmp(trialtimes_rift{r}.TRIAL_TYPE_RIFT,'RIFT1'))'
        if strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k},'negative_REFRAME') ...
                & strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k+1},'negative_REFRAME')
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k} = 'RtoR_1';
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k+1} = 'RtoR_2';
        elseif strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k},'negative_REFRAME') ...
                & strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k+1},'negative_AVOID')
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k} = 'RtoA_1';
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k+1} = 'RtoA_2';
        elseif strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k},'negative_AVOID') ...
                & strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k+1},'negative_REFRAME')
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k} = 'AtoR_1';
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k+1} = 'AtoR_2';
        elseif strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k},'negative_AVOID') ...
                & strcmp(trialtimes_rift{r}.ACTUAL_CONDITION{k+1},'negative_AVOID')
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k} = 'AtoA_1';
            trialtimes_rift{r}.TRIAL_TYPE_RIFT_PAIR{k+1} = 'AtoA_2';
        end
    end

    % Make names match ERT
    trialtimes_rift{r}.instructional_cue_started = trialtimes_rift{r}.instructional_cue_rift_started;
    trialtimes_rift{r}.instructional_cue_stopped = trialtimes_rift{r}.instructional_cue_rift_stopped;
    trialtimes_rift{r}.image_started = trialtimes_rift{r}.image_2_started;
    trialtimes_rift{r}.image_stopped = trialtimes_rift{r}.image_2_stopped;
    trialtimes_rift{r}.affect_rating_started = trialtimes_rift{r}.affect_rating_rift_started;
    trialtimes_rift{r}.affect_rating_stopped = trialtimes_rift{r}.affect_rating_rift_stopped;

    % Organize
    trialtimes_rift{r} = trialtimes_rift{r}(:, { ...
        'block_file', ...
        'image_category', ...
        'strategy', ...
        'new_strategy_chosen', ...
        'image_filename', ...
        'switched', ...
        'TRIAL_TYPE_ERT', ...
        'TRIAL_TYPE_RIFT', ...
        'TRIAL_TYPE_RIFT_PAIR', ...
        'ACTUAL_CONDITION', ...
        'instructional_cue_started', ...
        'instructional_cue_stopped', ...
        'image_started', ...
        'image_stopped', ...
        'affect_rating_started', ...
        'affect_rating_stopped' ...
        });

    writetable(trialtimes_rift{r},fullfile(inp.out_dir,['trialtimes_RIFT_run' num2str(r)]));

end


%% Find fmriprep files

% Scale motion params and save in SPM friendly format
for r = 1:2  % ert
    confD = dir([inp.(['fmriprep_ert' num2str(r) '_dir']) '/sub*/ses*/func/*_desc-confounds_timeseries.tsv']);
    conf = readtable(fullfile(confD(1).folder,confD(1).name),'FileType','text','Delimiter','tab');
    motT = conf(:,{'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'});
    mot = zscore(table2array(motT));
    writematrix(mot, fullfile(inp.out_dir,['motpar_ert' num2str(r) '.txt']))
end
for r = 1:4  % rift
    confD = dir([inp.(['fmriprep_rift' num2str(r) '_dir']) '/sub*/ses*/func/*_desc-confounds_timeseries.tsv']);
    conf = readtable(fullfile(confD(1).folder,confD(1).name),'FileType','text','Delimiter','tab');
    motT = conf(:,{'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'});
    mot = zscore(table2array(motT));
    writematrix(mot, fullfile(inp.out_dir,['motpar_rift' num2str(r) '.txt']))
end

% Find preprocessed image files, copy, and unzip
clear fmri_ert_nii
for r = 1:2
    niigzD = dir([inp.(['fmriprep_ert' num2str(r) '_dir']) '/sub*/ses*/func/*_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz']);
    fmri_ert_nii{r} = fullfile(inp.out_dir,['fmri_ert' num2str(r) '.nii']);
    copyfile( ...
        fullfile(niigzD(1).folder,niigzD(1).name), ...
        [fmri_ert_nii{r} '.gz'] ...
        );
end
clear fmri_rift_nii
for r = 1:4
    niigzD = dir([inp.(['fmriprep_rift' num2str(r) '_dir']) '/sub*/ses*/func/*_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz']);
    fmri_rift_nii{r} = fullfile(inp.out_dir,['fmri_rift' num2str(r) '.nii']);
    copyfile( ...
        fullfile(niigzD(1).folder,niigzD(1).name), ...
        [fmri_rift_nii{r} '.gz'] ...
        );
end
gunzip(fullfile(inp.out_dir,'fmri*.nii.gz'));
delete(fullfile(inp.out_dir,'fmri*.nii.gz'));

% Also the T1
niigzD = dir([inp.('fmriprep_ert1_dir') '/sub*/ses*/anat/*_space-MNI152NLin6Asym_desc-preproc_T1w.nii.gz']);
atlasT1_nii = fullfile(inp.out_dir,'t1.nii');
copyfile( ...
    fullfile(niigzD(1).folder,niigzD(1).name), ...
    [atlasT1_nii '.gz'] ...
    );
gunzip(fullfile(inp.out_dir,'t1.nii.gz'));
delete(fullfile(inp.out_dir,'t1.nii.gz'));


% Outputs for next step
outp = struct( ...
    'fmri_nii', {[fmri_ert_nii fmri_rift_nii]}, ...
    'trialtimes', {[trialtimes_ert trialtimes_rift]}, ...
    'hpf_sec', inp.hpf_sec, ...
    'fwhm_mm', inp.fwhm_mm, ...
    'out_dir', inp.out_dir ...
    );

