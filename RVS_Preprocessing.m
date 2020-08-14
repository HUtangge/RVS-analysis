%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   EEG analysis pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparing data
% .bdf and.sfp data for each subject

% 1. convert bdf. data into SPM format
% 2. Referencing (to common average)
% 3. High-pass filter
% 4. stop-filter (remove 50 Hz)
% 5. Downsamling to 512 Hz
% 6. Remove eye-blink 


%Experiment ID
file_mask = 'RVS';

% include some useful functions written by Jan Herding

% % add fieldtrip functions and use the default setting 
restoredefaultpath
addpath('D:\MATLAB\toolboxes\spm12');
addpath('D:\MATLAB\toolboxes\spm12\toolbox\MEEGtools');
addpath('D:\NoNobelPrizeWork\eeg_helper');
addpath('D:\MATLAB\toolboxes\spm12\external\fieldtrip');
ft_defaults

% define working directories
proj_dir    = 'D:\RVS_EEG';
rdata_dir   = fullfile(proj_dir, 'rawdata');
loc_dir     = fullfile(proj_dir, 'ElGuide');
%data_dir    = fullfile(proj_dir, 'preprocessed');

%sfp2mat(loc_dir, file_mask)
subj_ID = {'01', '02', '03', '04', '09', '11', '12', '13', '14', ...
          '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', ...
          '25' '26', '27', '28'};
      % exclusion : 10 due the high number of bad channels 
%% ========================================================================
%  1.) convert rawdata to SPM format & coregister electrode locations
% =========================================================================

for sj = 10:15
    clear D sj_dir rawfile outfiles
    %create a subject folder
    cd(proj_dir)
    sj_dir = fullfile(proj_dir, ['sj_' subj_ID{sj}]);
    if ~(exist(sj_dir, 'dir') == 7)
        mkdir(sj_dir);
    else
    end
    
    % find subject rawfile
    rawfile = fullfile(rdata_dir, [file_mask '_sj_' subj_ID{sj} '.bdf']);
    
    % name of the SPM file that'll be saved
    outfile = fullfile(sj_dir, file_mask);
    
    %read the EEG data & sensor location
    D = biosemi2spm(rawfile, outfile, [], loc_dir);
end

%%
for sj = 2:length(subj_ID)
    clear S D
    sj_dir = fullfile(proj_dir, ['sj_' subj_ID{sj}]);
    InterpData = fullfile(sj_dir, ['intp_RVS' subj_ID{sj} '.dat'])
    
    if exist(InterpData, 'file') == 2
        D = spm_eeg_load(InterpData);
        
        %% =============================================================
        %           2. RE-REFERENCE TO COMMON AVERAGE
        % ==============================================================
        % prefix of file name: M for montage
        S = [];
        S.D = D;
        S.refchan = 'average';
        D = spm_eeg_reref_eeg(S);
        
        %% =============================================================
        %           3. HIGH-PASS FILTER
        % ==============================================================
        % prefix of file name: f
        S = [];
        S.D = D;
        S.band = 'high';    % filterband [low|high|bandpass|stop]
        S.freq = 0.1;      % standard 0.01 - 0.1 for ERP, Cohen recommends 0.1-0.5
        %                    to the continuous data. He thinks bandpass filter is
        %                    most useful for time-frequency analysis
        S.type = 'butterworth';     %filter type [default: 'butterworth']
        %                           'butterworth': Butterworth IIR filter
        %                           'fir': FIR filter (using MATLAB fir1 function)
        D = spm_eeg_filter(S);
        
        %% ========================================================================
        %           4. Notch filter - removing 50 Hz noise
        % =========================================================================
        S = [];
        S.D = D;
        S.band = 'stop';
        S.freq = [49 51];
        S.type = 'butterworth';
        D = spm_eeg_filter(S);
    
        %% ========================================================================
        %                       5. Downsample to 512 Hz
        % =========================================================================
        % outputs files with prefix 'd'
        
        S = [];
        S.D = D;
        %choose a resampling method -   'resample' (default), 'decimate',
        %                               'downsample','fft'
        S.method = 'resample';
        S.fsample_new = 512;
        D = spm_eeg_downsample(S);
    else
    end
    clear sj_dir S D InterpData
end
%% ========================================================================
%                   6. remove Blink Artefacts                        
% =========================================================================
% make sure you have created a SPM EEG file that holds the spatial
% comppnent(s) reflecting a blink which you wanna remove from the EEG data
% outputs *ebf.*   eyeblink data

S = [];
S.D = D;

thresh = 11;    % change HERE until appropriate value for subject is found
                % approximately 100 blinks
D_ebf  = detect_eye_blinks(D, thresh); 

%% execute this part twice: Once to compute first two components as a check & finally only to compute the first component of eye-blink

% ==========================================================
%        compute (first) component of eye-blink
% ==========================================================
% fix the threshhold first and then use the D_ebf object with the correct threshold value to compute the eye-blink components
% first compute the first 2 principial components of the average eye-blink (num_components = 2) to check that only the first component captures reasonable
% electrodes involved in the eye-blink. If reasonable, recompute only the first principal component (num_components = 1), because we only want to remove the 
% first principal component later on. The results of compute_eye_blink_components is automatically saved as 'ebf_config' in the subject folder.

% outputs ebf_conf.mat and .dat file
num_components = 1;
compute_eye_blink_components(D_ebf, num_components);

%%
clear D_ebf
% clear old spatial confound structs
if isfield(D, 'sconfounds')
    S = [];
    S.D = D;
    S.method = 'clear';     %(EYES, BESA, SVD, SPMEEG, CLEAR)
    D = spm_eeg_spatial_confounds(S);
end

% add the spatial confound to the data set
S = [];clc
S.D = D;
S.method = 'SPMEEG';
S.conffile = 'ebf_conf.mat';
D = spm_eeg_spatial_confounds(S);

% remove eye blinks, outputs new files with prefix 'T'
S = [];
S.D = D;
S.correction = 'Berg';
D = spm_eeg_correct_sensor_data(S);
%% ========================================================================
%               7. Epoching
% ========================================================================
% determine the conditions
file_mask = 'RVS';
file_ext = 'StimLocked';
CurrPrefix = 'TdffMRVS';
proj_dir = 'D:\RVS_EEG';

subj_ID = '28';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the eye-blink corrected EEG file with prefix ^TdffM*
subj_dir = fullfile(proj_dir, ['sj_' subj_ID]);
cd(subj_dir)
continuous_data = fullfile(subj_dir, [CurrPrefix subj_ID '.dat']);
D = spm_eeg_load(continuous_data);

% copy file for epoched data in a new folder
EpochFolder = fullfile(subj_dir, file_ext);
if exist(EpochFolder, 'dir') ~= 7
    mkdir(EpochFolder)
else
end

cd(EpochFolder)

S         = [];
S.D       = D;
S.outfile = fullfile(EpochFolder ,[file_mask '_' file_ext '_' subj_ID]); % SPM 12
D = spm_eeg_copy(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find triggers that belong to the completed runs including in the analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = D.events;     %Onsets of the triggers


 % find out useful runs
    clear run_begin run_end run
    m = 1;
    n = 1;
    for j=1:(length(tmp)-1)
        % look for the first trigger after the event "Epoch", search
        % trigger 130 which suceeds an empty entry and is followed by
        % trigger 50 - 75
        if isempty(tmp(j).value) & ismember(tmp(j+1).value, 130) 
            if ismember(tmp(j+2).value, [50:75])
                run_begin(m) = j+1;
                m = m+1;
            end
        elseif (ismember(tmp(j).value, 130) & ismember(tmp(j+1).value, 130)) | (tmp(j).value > 9 & (ismember(tmp(j+1).value, 130)))
            % for the case in which we did not stop the recording between
            % two runs: look for two consecutive "130". 
            run_begin(m) = j+1;
            run_end(n) = j;
            m = m+1;
            n = n+1;
        elseif (tmp(j).value > 9 & isempty(tmp(j+1).value))
            if ~isempty(tmp(j-1).value)
                run_end(n) = j;
                n = n+1;
            end
        elseif (ismember(tmp(j).value, 130) & isempty(tmp(j+1).value))
            if ~isempty(tmp(j-1).value)
                run_end(n) = j;
                n = n+1;
            end
        end
        
    end
    if m ~= n
        run_end(end+1) = length(tmp);
    end
    
    % determine completed runs
    run = [run_begin' run_end'];
    run(:,3) = run(:,2) - run(:,1);
    run = run(run(:,3) > 490,:);   % exclude runs that have fewer than 490 entries
    run = run(valid_log{i}, :);


%%

% Replace empty enty with 999 (if the onsets of empty entries are equally
% spaced, they probably marks the start of a session 

for j=1:length(tmp)             
    if isempty(tmp(j).value)
        tmp(j).value=999;
    end
end

%get rid of the unvalid trial
%[tmp(1:20)] = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find the onsets of stimulation and rehearsal
% look for trigger 210, 220 230 and 240

for j = 1:length(tmp)
    % look for trigger 210, 220 230 and 240 and recode it to 333
    if tmp(j).value > 200 && mod(tmp(j).value, 10) == 0
       tmp(j).value = 333; 
    % look for trigger 10, 20, 30, 40, and recode to 555
    elseif tmp(j).value >= 10 && tmp(j).value <= 49 && mod(tmp(j).value, 10) == 0
         tmp(j).value = 555;
    end  
end

% IMPORTANT: HERE THE EVENTS ARE OVERWRITTEN BY THE RECODED TRIGGERS
D = events(D, 1,tmp);
save(D);

% define epochs
S    = [];
S.D  = D;
S.bc = 0; % baseline correction: off

epoch_interval = [-1000 5000];
S.timewin = epoch_interval;

% stimulation phase

S.trialdef.conditionlabel = 'Stimulation';
S.trialdef.eventtype      = 'STATUS';
S.trialdef.eventvalue     = 333;

S.reviewtrials = 0;
S.save         = 1;

%outputs epoched data 'eRVS_*'
D = spm_eeg_epochs(S);


%% ===========================================================
%  	           8. Automatic artefact removal
% ===========================================================

S                            = [];
S.D                          = D;
S.badchanthresh              = 0.2;
S.methods.channels           = 'EEG';
S.methods.fun                = 'jump';
S.methods.settings.threshold = 80; % µV
% S.methods.settings.twin      = artefact_twin;

D                            = spm_eeg_artefact(S);
%%           
S   = [];
S.D = D;
D   = spm_eeg_remove_bad_trials(S);
%%
% %% ===========================================================
% %  	           9. Low-pass filtering (only for ERP)
% % ============================================================
% 
S      = [];
S.D    = D;
S.band = 'low';
S.freq = 30;
D      = spm_eeg_filter(S);
% 
% 
% %% ===========================================================
% %  	           7. Baseline correction & averaging
% % ============================================================
% 
% 
% S         = [];
% S.D       = D;
% S.timewin = [-700 -200];
% D         = spm_eeg_bc(S);
% 
% S        = [];
% S.robust = 0;
% S.D      = D;
% S.review = 0;
% D = spm_eeg_average(S);
% % 
figure
hold on
color = {'r', 'm', 'g', 'b'}
for c = 1:4
plot(D.time, mean(D(D.indchannel({'CP4', 'CP6', 'P4', 'P6'}),:,c),1), 'Color', color{c}, 'LineWidth', 2)
end
% plot(D.time, mean(D(D.indchannel({'Fz', 'F2', 'FCz'}),:,1),1), 'k')
% 
% 
line([0 0], [-2 2], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
line([0.25 0.25], [-2 2], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
legend('f1', 'f2', 'f3', 'f4')
% line([1.25 1.25], [-3 3], 'Color', 'red', 'LineStyle', '--')
% line([2.50 2.50], [-3 3], 'Color', 'red', 'LineStyle', '--')
% line([3.75 3.75], [-3 3], 'Color', 'red', 'LineStyle', '--')
% 
% line([0.25 0.25], [-3 3], 'Color', 'red', 'LineStyle', '--')
% line([1.5 1.5], [-3 3], 'Color', 'red', 'LineStyle', '--')
% line([2.75 2.75], [-3 3], 'Color', 'red', 'LineStyle', '--')
% line([4.00 4.00], [-3 3], 'Color', 'red', 'LineStyle', '--')
% 
data=D.ftraw;
cfg=[];
data=ft_timelockanalysis(cfg,data)
% 
h=figure;
cfg=[];
ft_multiplotER(cfg,data)

%% Time-frequency transformation

%epoch_interval = [-1000 5000];
            %[4:2:48 52:4:100];    % in Hz [default: 4:2:48 Hz]
%TF_timewin     = (8./TF_frequencies) * 1000;       % only for multitaper (in ms) %%What's that??
%TF_timewin  = epoch_interval;
%TF_ncycles  = 7;        % Alex got the best result with 7, but try 5 first
%TF_timestep = 50;        % in ms, for a higher resolution 20ms 
%TF_freqres     = 0.2.*TF_frequencies; %[repmat(1./(TF_timewin./1000), 1, 24) ones(1,25).*5];         % in Hz only for multitaper
%TF_baseline    = [-700 -200];    % in sec f2-locked: [-2.25 -1.25]; saccades: whole interval [-1 1]; resp-locked: none ;) in ms for SPM12
%TF_prefix      = '2Hz'; % indicate finer resolution

S                   = [];
S.D                 = D;
S.channels          = {'EEG'};
S.frequencies       = 2:48;
S.timewin           = [-1000 5000];
S.phase             = 0;
S.method            = 'morlet';   %'mtmconvol'
S.settings.ncycles  = 7;
S.settings.subsample = 25; % round(temporal_resolution / (1000./D.fs)) in ms
%S.settings.timeres  = TF_timewin;
%S.settings.timestep = TF_timestep;
%S.settings.freqres  = TF_freqres;
%S.prefix            = TF_prefix;     
D                   = spm_eeg_tf(S);

%figure, imagesc(D.time, D.frequencies, squeeze(D(D.indchannel('F4'), :, :, 100)))
% axis xy
%%
D = spm_eeg_load('D:\RVS_EEG\sj_01\rehearsal_perf\tf_eRVS_reh_01.dat');
S                   = [];
S.D                 = D;
S.method            = 'LogR';
S.timewin           = [-250 -50];
S.Db                = 'D:\RVS_EEG\sj_01\baseline\tf_eRVS_base_01.dat';
D = ITIbaseline_spm_eeg_tf_rescale(S);
%%
h = figure(2);
for c = 1:4
    subplot(2,2,c)
    imagesc(D.time, D.frequencies, squeeze(mean(D(D.indchannel({'C4', 'C6', 'CP4', 'CP6'}),:,:,c),1)));
    axis xy
    xlim([-0.2 1.45])
    ylim([4 48])
    line([0 0], [4 48], 'Color', 'red', 'LineStyle', '--')
    line([0.25 0.25], [4 48], 'Color', 'red', 'LineStyle', '--')
    line([1.25 1.25], [4 48], 'Color', 'red', 'LineStyle', '--')
    colorbar
    hold on
end
%imagesc(D.time, D.frequencies, squeeze(mean(D(D.indchannel({'C3', 'C5', 'CP5', 'CP3','C4', 'C6', 'CP6', 'CP4', 'POz'}),:,:),1)));
axis xy
xlim([-0.2 1.45])
ylim([4 48])
line([0 0], [4 48], 'Color', 'red', 'LineStyle', '--')
line([0.25 0.25], [4 48], 'Color', 'red', 'LineStyle', '--')
line([1.25 1.25], [4 48], 'Color', 'red', 'LineStyle', '--')
% line([1.50 1.50], [4 48], 'Color', 'red', 'LineStyle', '--')    
% line([2.50 2.50], [4 48], 'Color', 'red', 'LineStyle', '--')
% line([2.75 2.75], [4 48], 'Color', 'red', 'LineStyle', '--')
% line([3.75 3.75], [4 48], 'Color', 'red', 'LineStyle', '--')
% line([4.00 4.00], [4 48], 'Color', 'red', 'LineStyle', '--')

colorbar



