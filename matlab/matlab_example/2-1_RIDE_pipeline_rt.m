%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 The RIDE Pipeline (matching with rt)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for instructions, see manual '#Allgemein\Manuale\EEG_Sprachartefakt_Korrektur\manual\Manual_RIDE_pipeline.docx' 

% This script is used for the preparation of single-trial EEG data, after RIDE correction (pre-processing in BVA) 
% the EEG data will be sorted in a specific order, which can be used to concatenate it with the log file 
% it will generate a folder named 'export_R' for later analysis (e.g., LMM in R)

% data preparation
% Before you start: this script will access 9 different subfolders in your experiments RIDE folder: 
%"raw_eeg", "raw", "mat", "BVAinfo", "output", "single_trials", "single_trials_small", "log", "log_c", "final_log" and "export_R". Please create them now. 

%NEW folders:
% raw_eeg: raw data before BVA has to be saved here; will be used to extract the trigger information if RT is used to match the data (EEG + log).
% raw: MAT-files created by the export function in BVA have to be saved here (see former RIDE script). Matlab does not write here. 
% mat: matlab created (imported) data files (during RIDE), including e.g. reaction times % read from your voice triggers. 
% output: matlab created: results from the actual RIDE calculations. 
% BVAinfo: matlab created: the information from EEG data after BVA, which can be used for check with log file.
% single_trials: matlab created: single trial EEG data after RIDE correction, and the data with interested time window.
% single_trials_small: make the data set samll
% log: this folder is the original log files
% log_c: matlab created: this folder is the log files with .mat format, it only reads numerial data from log file
% final_log: matlab created: this has the log information matched with eeg data
% export_R: matlab created: this has the information of interested timewindow and the text file generated from these files

%

% Tang Ge


clear all;
% add toolbox used in this script
addpath(genpath('N:\Software\Matlab_toolboxen\eeglab13_5_4b'));
addpath(genpath('N:\matlab\'));%'N:\Software\Matlab_toolboxen\RIDE_call'));

% add working folder for this project
addpath(genpath('C:\Users\neuro-lab\Documents\Research projects\eeg-ride\matlab\matlab_example\'));

%% Basic Path of your experiments RIDE folder (this folder contains your raw, mat and output folders): former RIDE script (Guang)

% define folders for the processed data 
ridefolder = 'C:\Users\neuro-lab\Documents\Research projects\eeg-ride\matlab\matlab_example\';

twd = [-100,1200]; %the time window for the epoched data
        % NOTE: if segmenated epoch is too long and covers triggers from the next trial, make sure
        % to run the 'segmentation' in BVA again to make it ONLY covers two triggers! [Pei, May2020]
         
sub = {'Vp0001'};

% Your condition names - Bedingungsnamen wie sie im filename enthalten sind
con = {'rep1_close','rep1_distant','rep1_het',...
       'rep2_close','rep2_distant','rep2_het', ...
       'rep3_close','rep3_distant','rep3_het', ...
       'rep4_close','rep4_distant','rep4_het', ...
       'rep5_close','rep5_distant','rep5_het'};
        % NOTE: it is crucial that the sum of these data cover *exactly*
        % the trial numbers of your experiment. See "1-1_Manual". [Pei, July 2020]

% Voicekeytrigger 
voice_tr = 'S  4';

% condition trigger, which has been used for segmentation in BVA
contr = {'S201','S202','S203'};

% sampling rate
srate = 500;


%%
% %% Preparation for RIDE correction and matching;
% % ---------------------------------------
% % ONLY NEEDS TO BE DONE ONCE IN THE BEGINNING
% % get the EEG data and reaction time (needed for RIDE), which is being extracted from EEG data (former RIDE script, Guang)
% 
% % create new folders or, if they have already beeen created, check if they are there
% mat_dir = fullfile(ridefolder, 'mat');
% if exist(mat_dir,'dir')~=7
%     mkdir(mat_dir)
% end
% BVA_dir = fullfile(ridefolder, 'BVAinfo');
% if exist(BVA_dir,'dir')~=7
%     mkdir(BVA_dir)
% end
% 
% for j = 1:length(sub) %loops through subjects
%     disp(j); % display subject id
%     
%     for k = 1:length(con) % loops through condition names to later read files in
%         
%         % this is load data from export file after BVA
%         EEG = pop_loadbva([ridefolder, 'raw\', sub{j},'_',con{k},'.mat']);
%         EEG.setname=[sub{j},'_',con{k}]; % set setname using subject id and condition
%         EEG = eeg_checkset( EEG ); % look at EEG
%         
%         rt = nan(size(EEG.data,3),1); %creates new array for RTs 
%         triggers = {EEG.event.type}; % read out all triggers in data set;  TEST: EEG.event(3).type
%         epochs = [EEG.event.epoch];
%         latency = [EEG.event.latency];
%         contr_value_all = []; % create a new variable in which we insert unique trigger (created above) from BVA (EEG segments)
%         
%         for t = 1:size(EEG.data,3) % loops through all trials (per condition)
%             
%             temp = triggers(epochs==t); % extract all triggers of one EEG segment
%             temp1 = latency(epochs==t); % stores trigger position (in ms); needed to extract RT latency later
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            
%             %sometimes, there are multiple condition or voice triggers  in
%             %one segment due to overlapping segments;
%             %here, we only extract the FIRST condition and the FIRST voice
%             %trigger of each segment
%            
%             if any(strcmpi(temp,voice_tr)) && any(ismember(temp,contr)) % if there is a voice trigger and a condition trigger MATlab continues processing that segment
%                 
%                 voice_tr_ind = find(strcmpi(temp,voice_tr));
%                 tl = temp1(voice_tr_ind(find(strcmpi(temp,voice_tr)) > find(ismember(temp,contr),1))); %in case of multiple response triggers (voice key triggert multiple times) get the first trigger and then next line read out only first response latency
%                 cts = temp1(ismember(temp,contr));
%                 
%                 % Tang Ge: to get the triggers
%                 contr_value = temp(ismember(temp, contr));
%                 contr_value = cell2mat(cellstr(contr_value));  %% defined above
%                 contr_value = regexp(contr_value, '\d*','match'); 
%                 contr_value = str2double(contr_value);
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % RTs are being calculated
%                 % optional: remove trials with too short RTs 
%                 
%                 if ~isempty(tl)
%                     % tl(1) = first response trigger; cts(1) = first condition trigger; change numbers if necessary
%                     rt(t) = tl(1) - cts(1);% read latencies out for trials containg both condition AND response triggers and subtract from each other to have CORRECT (nor recording) latencies
%                 else
%                     rt(t) = nan(1);
%                 end
%                 
% %                 % remove trials with RTs less than 500 ms (i.e. 250 data points)
% %                 if rt(t) < 250  % optional
% %                     rt(t) = nan(1); % optional
% %                 end   %optional
%                 
%             end
%             contr_value_all = [contr_value_all; contr_value];
%             
%         end
%                 
%         contr_value_all = contr_value_all(~isnan(rt),:); % gets the condition triggers
%         data = EEG.data(:,:,~isnan(rt)); % gets "real" data read out from EEG.data, but only those trials that don't have missing (or too short) RTs
%         rt = rt(~isnan(rt)); % now make new RT vector containing only available RTs (remove still "empty" [=NAN] cells)
%         data = permute(data,[2,1,3]); % changes structure: electode * time * trial -> time * electrode * trial
%         rt = rt*1000/EEG.srate; % coverts to unit of millisecond
%         save([ridefolder,'mat/',sub{j},'_',con{k},'.mat'],'data','rt'); % mat_control= new folder ;;save file % Path definitions correct?%%'mat\',sub{j},'_',con{k},'.mat'],'data','rt');
%         
%         %TangGe: infoBVA : info about RT and trigger to merge files later
%         infoBVA = [rt contr_value_all];
%         save([ridefolder,'BVAinfo/', sub{j},'_',con{k},'_infoBVA.mat'], 'infoBVA'); % saves the information we get from EEG data after BVA and put them in folder infoBVA
%         
%     end
% end
% 
% chanlocs = EEG.chanlocs; % saves channel locations needed for later plotting
% save([ridefolder,'chanlocs.mat'],'chanlocs');
% % needs to be done only once, otherwise overwritten every time
% 
% % clear useless variables
% clear mat_dir BVA_dir % j k t temp temp1 tl cts

%% RIDE - One subject, one condition
output_dir = fullfile(ridefolder, 'output');
j = 1;
k = 1;
  load([ridefolder,'mat\',sub{j},'_',con{k},'.mat'],'data','rt');    % Path definitions correct?
  
        cfg = [];
     % Abtastrate ist 500Hz, daher 2, da alle 2ms aufgenommen wird
        cfg.samp_interval = 2;
     % Segmentlänge (geht z.B. bis 2000ms aber der letzte Punkt der aufgenommen wird ist 1998 (2000 - 2)
        cfg.epoch_twd = [-100,1198];
     % nicht ändern
        cfg.comp.name = {'s','r'}; %'c': based on visual inspection (Usually a component is identified by a big and wide hump in the late time window )
     % response trigger
        cfg.comp.twd = {[0,600],[-300,300]}; % rt set according to RIDE_implementtation.pptx p.15
            % time window parameters (for s und r component) have to be adapted to experiment, 
            % e.g., PWI, compounds (ANTJE): [0, 550], [-250,1000] / AK: 0:800; -600:1000
            % BSD: (17 June 2020) [0,800],[-600,1000]
     % cfg.comp.latency = {0,'unknown'}; for c component
        cfg.comp.latency = {zeros(size(data,3),1),rt}; % BSD:(17 June 2020) {0,rt}
        cfg = RIDE_cfg(cfg);
        results = RIDE_call(data,cfg);
        save([ridefolder,'output\results_',sub{j},'_',con{k},'.mat'],'results'); % output folder


% %% RIDE calculations Guang script
% % check all parameters in this section and change if necessary.
% 
% % 21. May 2020: error: probably due to the unremoved A1 channel
% % -> rerun BVA analysis and generate new workspace for the 62-channels
% % 22. May 2020: error: last step errornously compute 2-D data
% % -> rerun BVA analysis and re-segment to exclude overlapping triggers
% 
% output_dir = fullfile(ridefolder, 'output');
% if exist(output_dir,'dir')~=7
%     mkdir(output_dir)
% end
% 
% for j = 1:length(sub)
%     disp(j);
%     
%     for k = 1:length(con)
%         
%         load([ridefolder,'mat\',sub{j},'_',con{k},'.mat'],'data','rt');    % Path definitions correct?
%         cfg = [];
%      % Abtastrate ist 500Hz, daher 2, da alle 2ms aufgenommen wird
%         cfg.samp_interval = 2;
%      % Segmentlänge (geht z.B. bis 2000ms aber der letzte Punkt der aufgenommen wird ist 1998 (2000 - 2)
%         cfg.epoch_twd = [-100,1198];
%      % nicht ändern
%         cfg.comp.name = {'s','r'}; %'c': based on visual inspection (Usually a component is identified by a big and wide hump in the late time window )
%      % response trigger
%         cfg.comp.twd = {[0,600],[-300,300]}; % rt set according to RIDE_implementtation.pptx p.15
%             % time window parameters (for s und r component) have to be adapted to experiment, 
%             % e.g., PWI, compounds (ANTJE): [0, 550], [-250,1000] / AK: 0:800; -600:1000
%             % BSD: (17 June 2020) [0,800],[-600,1000]
%      % cfg.comp.latency = {0,'unknown'}; for c component
%         cfg.comp.latency = {zeros(size(data,3),1),rt}; % BSD:(17 June 2020) {0,rt}
%         cfg = RIDE_cfg(cfg);
%         results = RIDE_call(data,cfg);
%         save([ridefolder,'output\results_',sub{j},'_',con{k},'.mat'],'results'); % output folder
%     
%     end
% end
% 
% clear output_dir j k 

%% simply plot the separation senario

for j = 1:length(sub)
    for k = 1:length(con)
        load([ridefolder,'output\results_',sub{j},'_',con{k},'.mat'],'results'); % Attention! Folder here must match folder in line 109
        erp_all(:,:,j,k) = results.erp;
        erp_new_all(:,:,j,k) = results.erp_new;
        s_all(:,:,j,k) = results.s;
%         c1_all(:,:,j,k) = results.c1; 
%         c2_all(:,:,j,k) = results.c2; 
        r_all(:,:,j,k) = results.r; % speech artifact
    end
end

clear j k 

t_axis = linspace(twd(1),twd(2), (twd(2) - twd(1))*srate/1000); % old length 1400 --> interval and x-line must agree

figure;

subplot(1,3,1);
plot(t_axis,mean(mean(erp_all,4),3));
axis([twd(1),twd(2),-20,20]);
xlabel('time');
title('erp');
hold on;

% subplot(1,3,2);
% plot(t_axis,mean(mean(erp_new_all,4),3));
% axis([twd(1),twd(2),-20,20]);
% xlabel('time');
% title('erp_new');
% hold on;

subplot(1,3,2);
plot(t_axis,mean(mean(s_all,4),3));
axis([twd(1),twd(2),-20,20]);
xlabel('time');
title('s');
hold on;

subplot(1,3,3);
plot(t_axis,mean(mean(r_all,4),3));
axis([twd(1),twd(2),-20,20]);
xlabel('time');
title('r');

% subplot(1,3,2);
% plot(t_axis,mean(mean(c1_all,4),3));
% axis([twd(1),twd(2),-20,20]);
% xlabel('time');
% title('c1');
% hold on;
% 
% subplot(1,3,3);
% plot(t_axis,mean(mean(c2_all,4),3));
% axis([twd(1),twd(2),-20,20]);
% xlabel('time');
% title('c2');
% hold on;

%% Useful plottings (From Guang's RIDE manual)
chan_index = find(strcmpi({chanlocs.labels},'Pz')); % get channel infos

% Example 1: show the ERP and the RIDE components for this subject:
figure;RIDE_plot(results,{'erp','s','c1','c2','r'},chan_index); 

% Example 2: show the ERP the re-constructed ERP:
figure;RIDE_plot(results,{'erp','erp_new'},chan_index);

% Example 3: show the topographies (needs EEGLAB)
t = 500;t1 = fix((t-cfg.epoch_twd(1))/cfg.samp_interval);
figure;subplot(1,5,1);topoplot(results.erp(t1,:),chanlocs);text(0,1,'erp');
subplot(1,5,2);topoplot(results.s(t1,:),chanlocs);text(0,1,'s');
subplot(1,5,3);topoplot(results.c1(t1,:),chanlocs);text(0,1,'c1');
subplot(1,5,4);topoplot(results.c2(t1,:),chanlocs);text(0,1,'c2');
subplot(1,5,5);topoplot(results.r(t1,:),chanlocs);text(0,1,'r');

%% Guang script, adapted to be used for all participants in parallel by Tang Ge: 
% script extracts single trial data from (aggregated) Ride corrected data;
% 3-dimensional matrix file + 2-dimensional results file are created during RIDE correction (s. mat and output folder)
% this script: R - component(per participant & condition) is subtracted from
% each single trial (mat files: 3-d matrix before RIDE-correction --> cleaned single trials for statistical analysis;
% this part is for one subject; 
% create new folder: single_trials 
%%%%%%%% If you have already run RIDE correction and do not use additional
%%%%%%%% trigger information to match your files then you start your
%%%%%%%% analyses here. (AK, April 2019)

single_trials_dir = fullfile(ridefolder, 'single_trials');
if exist(single_trials_dir,'dir')~=7
    mkdir(single_trials_dir)
end

for j = 1:length(sub) % loop over subjects
    disp(j);
    single_trial_S_all = [];
    single_trial_R_all = [];
    single_trial_erp_all = []; % erp before RIDE
%     single_trial_erp_new_all = []; % re-constructed erp
    
    for k = 1:length(con) % loop over conditions
        load([ridefolder,'mat/',sub{j},'_',con{k},'.mat']); %3-dimensional matrix, data before Ride correction
        load([ridefolder,'output/results_',sub{j},'_',con{k},'.mat']); % 2-dimensional file, after Ride correction
        
        % get single trial S (free of articulation artifact)
        single_trial_S = data - move3(results.r(:,:,ones(length(rt),1)),round(results.latency_r*srate/1000)); % remove R
%         single_trial_R = data - results.s(:,:,ones(length(rt),1)); % remove S
        single_trial_erp = data; % original erp
%         single_trial_erp_new = data - results.c1(:,:,ones(length(rt),1)) - results.c2(:,:,ones(length(rt),1))...
                                - move3(results.r(:,:,ones(length(rt),1)),round(results.latency_r*srate/1000)); % remove Cs & R
%         % alternatively: c1, c2 also latency-locked
%         single_trial_erp_new = data - move3(results.c1(:,:,ones(length(rt),1)),round(results.latency_r*srate/1000))...
%             - move3(results.c2(:,:,ones(length(rt),1)),round(results.latency_r*srate/1000))...
%             - move3(results.r(:,:,ones(length(rt),1)),round(results.latency_r*srate/1000)); 
                            
        % reshape the data for importing to 'export function'
        single_trial_S = permute(single_trial_S,[2,1,3]);
%         single_trial_R = permute(single_trial_R,[2,1,3]);
        single_trial_erp = permute(single_trial_erp,[2,1,3]);
%         single_trial_erp_new = permute(single_trial_erp_new,[2,1,3]);
        
        single_trial_S_all = cat(3, single_trial_S_all, single_trial_S);
%         single_trial_R_all = cat(3, single_trial_R_all, single_trial_R);
        single_trial_erp_all = cat(3, single_trial_erp_all, single_trial_erp);
%         single_trial_erp_new_all = cat(3, single_trial_erp_new_all, single_trial_erp_new);
        
    end
    
    % save all the single trial data to single_trials, subject by subject
    save([ridefolder,'single_trials/',sub{j},'_S.mat'], 'single_trial_S_all'); 
%     save([ridefolder,'single_trials/',sub{j},'_R.mat'], 'single_trial_R_all');
    save([ridefolder,'single_trials/',sub{j},'_erp.mat'], 'single_trial_erp_all');
%     save([ridefolder,'single_trials/',sub{j},'_erp_new.mat'], 'single_trial_erp_new_all');
    
    clear single_trial_S single_trial_erp single_trial_S_all single_trial_erp_all 
%     clear single_trial_R single_trial_erp_new single_trial_R_all single_trial_erp_new_all 
    
end

clear single_trials_dir j k 

%% making your dataset smaller by using a reduced time window 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put the reduced time window here
new_twd = [-100 600]; % v1:900
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
small_dir = fullfile(ridefolder, 'single_trials_small');    % ?
if exist(small_dir,'dir')~=7
    mkdir(small_dir)
end

times_all = twd(1):(1000/srate):twd(2); % original timewindow
times_selected = new_twd(1):(1000/srate):new_twd(2); % time window small
[times_shared, times_idx] = intersect(times_all, times_selected, 'stable');
times_idx = times_idx';

for j = 1:length(sub)
    
    load([ridefolder,'single_trials/',sub{j},'_S.mat']); % load the data
%     load([ridefolder,'single_trials/',sub{j},'_R.mat']);
    load([ridefolder,'single_trials/',sub{j},'_erp.mat']);
%     load([ridefolder,'single_trials/',sub{j},'_erp_new.mat']);
    
    single_trial_S_times_selected = single_trial_S_all(:, times_idx, :);
%     single_trial_R_times_selected = single_trial_R_all(:, times_idx, :);
    single_trial_erp_times_selected = single_trial_erp_all(:, times_idx, :);    
%     single_trial_erp_new_times_selected = single_trial_erp_new_all(:, times_idx, :);
    
    save([ridefolder,'single_trials_small/',sub{j},'_S_small.mat'], 'single_trial_S_times_selected'); % save all the single trial data to single_trials, subject by subject
%     save([ridefolder,'single_trials_small/',sub{j},'_R_small.mat'], 'single_trial_R_times_selected');
    save([ridefolder,'single_trials_small/',sub{j},'_erp_small.mat'], 'single_trial_erp_times_selected');
%     save([ridefolder,'single_trials_small/',sub{j},'_erp_new_small.mat'], 'single_trial_erp_new_times_selected');
    
    delete([ridefolder,'single_trials/',sub{j},'_S.mat']); % deletes large data file to save space 
%     delete([ridefolder,'single_trials/',sub{j},'_R.mat']); 
    delete([ridefolder,'single_trials/',sub{j},'_erp.mat']); 
%     delete([ridefolder,'single_trials/',sub{j},'_erp_new.mat']); 
end

clear small_dir j single_trial_S_times_selected single_trial_erp_times_selected 
% clear single_trial_R_times_selected single_trial_erp_new_times_selected

%% Generate one file for all subject, which can be used for later ERP analysis

S_all = []; % stimulus-locked (no R)
% R_all = []; % response-locked
erp_all = []; % BVA preprocessed, segmented data before RIDE
% erp_new_all = []; % stimulus-locked (no C1,C1,R; pure S)

for j = 1:length(sub)
    
    % if you use the small data set, then you need to change to the folder with samll data set
    load([ridefolder,'single_trials_small/',sub{j},'_S_small','.mat']);
%     load([ridefolder,'single_trials_small/',sub{j},'_R_small','.mat']);
    load([ridefolder,'single_trials_small/',sub{j},'_erp_small','.mat']);
%     load([ridefolder,'single_trials_small/',sub{j},'_erp_new_small','.mat']);
    
    S_all = cat(3, S_all, single_trial_S_times_selected);
%     R_all = cat(3, R_all, single_trial_R_times_selected);
    erp_all = cat(3, erp_all, single_trial_erp_times_selected);
%     erp_new_all = cat(3, erp_new_all, single_trial_erp_new_times_selected);
end
save([ridefolder,'single_trials_small/S_all.mat'], 'S_all', '-v7.3');% save all the EEG data for all participants and all conditions in one big file 
% save([ridefolder,'single_trials_small/R_all.mat'], 'R_all', '-v7.3');
save([ridefolder,'single_trials_small/erp_all.mat'], 'erp_all', '-v7.3');
% save([ridefolder,'single_trials_small/erp_new_all.mat'], 'erp_new_all', '-v7.3');

clear j

%% Here, log files are being prepared for matching                    
% Read txt file (log files of each participant) into mat format
% this part only reads numbers into mat file (per column)

log_dir = fullfile(ridefolder, 'log_c');   
if exist(log_dir,'dir')~=7
    mkdir(log_dir)
end

for i = 1: length(sub) %1: length(sub)
    ind = str2double(cell2mat(regexp(sub{i}, '\d*','match'))); % get the index of subject    
        % Note: rename raw log files, e.g. Vp0001.txt -> 1.txt [Pei, May 2020]
    fileID = fopen([num2str(ind) '.txt'], 'r'); % filename = fopen(fileID) --> check file path
    log = [];
    line = fgetl(fileID);  % this is the headline, we do not need to put in log file
    line = fgetl(fileID); % this is the empty line, we do not need to put in log file [Pei, May 2020]
    line = fgetl(fileID); % this is the official first line of data [Pei, May 2020]
    while 1
        log_single = regexp(line,'\s','split'); % get all the numerical data from txt file  
        log_single = log_single(~cellfun('isempty',log_single));
        log_single = log_single([1 4 9 10 11]); % refers to columns in log-file that we want to extract; 
            % in this example the column numbers for: 
            % "subject","item_ID", "RT", "blocking condition", "repetition cycle"
            % Note: make sure subject variable contains ONLY number! [Pei, May 2020]
        log_single = str2double(log_single);
        log_single = log_single(~isnan(log_single));
        log = [log; log_single];
        line = fgetl(fileID); % move to next line       
        if line == -1
            break
        end
    end
    fclose(fileID);
    log(:,end+1) = [1:length(log)]; % give each trial an index
    save([ridefolder, 'log_c\' sub{i}, '_log.mat'], 'log'); % save to specific direction; log_c is a new folder, which had been created before
end

clear log_dir i ind line fileID log_single log

%% Read triggers from raw eeg
% #BSD ==> note overlapping trigger names! [Pei, June2020]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an example from Pei's #BSD experiment
% trigger settings for the experiment, it should have all the triggers used in the experiment
trigger.item_id = [1:125];
trigger.repetition = [1:5];
% this trigger is used to define segment(epoch) it is better to keep this name: trigger.con_trigger
trigger.con_trigger = [201 202 203]; 
% this trigger is also mandatory, if it can be voice trigger or response trigger.
trigger.response = [4]; % in this Anna K's case, the response trigger is the voice trigger
trigger.error = [254];

% this is important, the position of condition trigger
con_trigger_pos = 3; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the variable for the import of eeg_readings
% check if eeg_log folder already exists or create one
eeg_log_dir = fullfile(ridefolder, 'eeg_log');
if exist(eeg_log_dir,'dir')~=7
    mkdir(eeg_log_dir)
end

fn = fieldnames(trigger);

for sj = 1:length(sub)
    
    % read trigger information from raw eeg file with fieldtrip
    headerfile = [sub{sj} '.vhdr'];
    EEG = pop_loadbv(fullfile(ridefolder, 'raw_eeg'), headerfile);
    
    % get epoch information trial by trial
    tmp = {EEG.event.type}'; % read all triggers
    for i = 1:length(tmp)
        if isempty(tmp{i}) | ismember(tmp{i}, 'boundary') % for missing & boundary input
            tmp{i} = 'S 9999'; % fill in all the triggers
        end
    end
    clear i
    tmp = cellstr(tmp);
    tmp = regexp(tmp, '\d*', 'match'); % 1st column: trigger (only number info)
    tmp = str2double([tmp{:}])'; 
    tmp(:,2) = [EEG.event.latency]'; % 2nd: latency in raw eeg
    tmp(:,3) = [1:size(tmp,1)]; % 3nd: index in raw eeg
    tmp(find(tmp(:,1) == 9999), :) = [];
    
    % create an empty evtlog matric
    con_triggers = tmp(ismember(tmp(:,1), trigger.con_trigger),1);
    evtlog = zeros(length(con_triggers), (length(fn) + 2)); 
             % the first fiive columns: trigger infos (item, rep, blocking, response, errorkey)
             % the last two columns are reaction time and the position of these triggers in the raw eeg data
    evtlog(:, con_trigger_pos) = con_triggers; % assign con_trigger to the contrast trigger positioin (3nd)
    tmp(end+10, :) = 0; % add 10 empty rows in the end (?)
  
    j = 0;
    for i = 1:length(tmp) % read in each line in raw eeg
        if ismember(tmp(i,1), trigger.con_trigger) % if the trigger refers to a segmentation info (blocking)...
            j = j+1; % ... add new trial
           
            if i-con_trigger_pos >= 0 % except the very first two lines
                for k = 1 % column 1: itemID
                    evt_ind = find(ismember(tmp((i-con_trigger_pos+1),1), trigger.(fn{k}))) + i - con_trigger_pos;
                    % finds in the raw eeg whether "2 rows above the con_trigger row" (i.e. trigger item ID)
                    % matches the item trigger, and assign in the new eeg_log an index
                    if ~isempty(evt_ind)
                        evtlog(j,k) = tmp(evt_ind,1); % trial j, column 1
                    end
                end
                for k = 2 % column 2: repetition
                    evt_ind = find(ismember(tmp((i-1),1), trigger.(fn{k}))) + i - con_trigger_pos + 1;
                    % finds in the raw eeg whether "1 rows above the con_trigger row" (i.e. trigger repetition)
                    % matches the repetition trigger, and assign in the new eeg_log an index
                    if ~isempty(evt_ind)
                        evtlog(j,k) = tmp(evt_ind,1); % trial j, column 2
                    end
                end
                for k = 4 % column 4: response
                    evt_ind = find(ismember(tmp((i+1),1), trigger.(fn{k}))) + i;
                    % finds in the raw eeg whether "1 rows below the con_trigger row" (i.e. trigger response)
                    % matches the response trigger, and assign in the new eeg_log an index
                    if ~isempty(evt_ind)
                        evtlog(j,k) = tmp(evt_ind,1); % trial j, column 4
                        if ismember(evtlog(j,k), trigger.response) % if the trigger refers to a response
                            evtlog(j,length(fn)+1) = tmp(evt_ind,2) - tmp(i,2);
                            % column 6: calculate RT by substracting time of voice key from time of con_trigger
                        end
                    end
                end
                for k = 5 % column 5: error key
                    evt_ind = find(ismember(tmp((i + length(fn) - con_trigger_pos),1), trigger.(fn{k}))) + i + 1;
                    % finds in the raw eeg whether "2 rows below the con_trigger row" (i.e. trigger error key)
                    % matches the  error key trigger, and assign in the new eeg_log an index
                    if ~isempty(evt_ind) % if correct (not marked as error)
                        evtlog(j,k) = tmp(evt_ind,1); % trial j, column 5
                    end
                end
                
            else % for the very first two lines
                for k = 1:(con_trigger_pos-1)
                    evt_ind = find(ismember(tmp(1:(i-1),1), trigger.(fn{k})), 1, 'first');
                    if ~isempty(evt_ind)
                        evtlog(j,k) = tmp(evt_ind,1);
                        if evtlog(j,k) == trigger.response
                            evtlog(j,length(fn)+1) = tmp(evt_ind,2) - tmp(i,2);
                        end
                    end
                end
                for k = (con_trigger_pos+1):length(fn)
                    evt_ind = find(ismember(tmp((i+1): (i + length(fn) - con_trigger_pos),1), trigger.(fn{k})), 1, 'first') + i;
                    if ~isempty(evt_ind)
                        evtlog(j,k) = tmp(evt_ind,1);
                        if evtlog(j,k) == trigger.response
                            evtlog(j,length(fn)+1) = tmp(evt_ind,2) - tmp(i,2);
                        end
                    end
                end
            end
            evtlog(j,length(fn)+2) = tmp(i,3); % ...column 7: position of the triggers in the raw eeg data
        end
        
    end

    clear i j k
    
    save(fullfile(ridefolder, 'eeg_log', [sub{sj} '_eeg_log.mat']), 'evtlog'); % save the file in folder eeg_log
    clear evtlog tmp
    
end

%% Connect information from raw eeg file with behavioral file

% these are triggers used for matching eeg_log and log_c, you can pick triggers have the most information in your experiment
% if you want to add more triggers, you can do with line 475
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the condition trigger, trigger used for segment or epoch is mandatory
evt_con_trigger_pos = 3; % the trigger position in EEG file
log_con_trigger_pos = 4; % the trigger position in log file
evt_picid_pos = 1;
log_picid_pos = 2;
evt_rep_pos = 2;
log_rep_pos = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = [1:16 18:24] % 1:length(sub) % (1)17; (2)[1:16 18:24]
    % Note: Vp0018 only has 1275 trials in EEG recordings, run separately
    % [Pei, June 2020 #BSD]
    load([ridefolder, 'eeg_log\', sub{i}, '_eeg_log.mat']);
    load([ridefolder, 'log_c\', sub{i}, '_log.mat']);
    ind = [1:size(evtlog,1)];
   
    if length(ind) == size(log,1)
        if sum(evtlog(:,evt_con_trigger_pos) - log(:,log_con_trigger_pos)) == 0
            log = [log evtlog ind'];
        else
            error('Missing trials');
        end
    else
        for k = 1:(size(log,1)-2)          
            % this line you can have more conditions to get a precise match
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if k <= size(evtlog,1) & ...
                    ((evtlog(k,evt_picid_pos)~=0 & evtlog(k,evt_picid_pos) ~= log(k,log_picid_pos)) & ...% this line does not match
                   (evtlog(k+1,evt_picid_pos) == log(k+1,log_picid_pos))) ... % but if the next line matches
                |  ((evtlog(k,evt_picid_pos)~=0 & evtlog(k,evt_picid_pos) ~= log(k,log_picid_pos)) & ...
                   (evtlog(k+2,evt_picid_pos) == log(k+2,log_picid_pos)))     % but if the next next line matches
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                evtlog(k,:) = nan(1); % mark ass NAN
            elseif k <= size(evtlog,1) & ...
                    (evtlog(k,evt_picid_pos)~=0 & evtlog(k,evt_picid_pos) ~= log(k,log_picid_pos)) % this line does not match
                evtlog((k+1):(end+1),:) = evtlog(k:end,:); % insert new line
                evtlog(k,:) = nan(1); % mark as NAN (missing data)
            end
        end
        
       for k = (size(log,1)-1):size(log,1) % the last rows    
           if k <= size(evtlog,1) & ...
                    (evtlog(k,evt_picid_pos)~=0 & evtlog(k,evt_picid_pos) ~= log(k,log_picid_pos)) % this line does not match
                evtlog((k+1):(end+1),:) = evtlog(k:end,:); % insert new line
                evtlog(k,:) = nan(1); % mark as NAN (missing data)
            end
        end
        
        if size(evtlog,1) < size(log,1)
            last_trial = size(evtlog,1);
            evtlog((end+1):size(log,1), :) = nan(1); % add NANs
        elseif size(evtlog,1) > size(log,1)
            evtlog(size(log,1)+1:end, :) = []; % delete final rows
        end
         
        ind = [1:size(evtlog,1)];
        log = [log evtlog ind'];
        
    end
    
    save([ridefolder, 'log_c\' sub{i}, '_log.mat'], 'log');
    clear log evtlog ind
    
end

clear i k

%% MATCHING EEG and Behavioral data 

final_log_dir = fullfile(ridefolder, 'final_log');
if exist(final_log_dir,'dir')~=7
    mkdir(final_log_dir)
end

for i = 1:length(sub) %1:length(sub) % note: 16 = vp17 
    load([ridefolder, 'log_c\', sub{i}, '_log.mat']);
    final_log = [];
    
    for j = 1:length(con) % 1:length(con)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This part need to be change according to different experiemnts 
        % find conditions in the log file (cause they are not sorted)
        
        % repetition cycle
        rep = str2double(regexp(con{j}, '\d', 'match')); % rep 1,2,3,4,5
        if isempty(rep) % all cycles
            rep = [1:5];
        end

        % blocking condition
        if ismember('c', con{j}) 
            blocking = 201; % close
        elseif ismember('d', con{j}) 
            blocking = 202; % distant
        elseif ismember('h', con{j}) 
            blocking = 203; % het
        else 
            blocking = [201:203]; % all conditions
        end
        
        % this code should catch all the trials in one condition
        cond_ind = find(ismember(log(:,4), blocking) & ismember(log(:,5), rep)); % index column position in log
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % load the file from BVA_info
        load([sub{i},'_',con{j},'_infoBVA.mat']);
        
        % select the information for matching from log file
        log_ind = log(cond_ind, end);  % extract the trial index for each condition, you should find which column is trial index
        log_rt = log(cond_ind,end-2);  % extract the special trigger for each condition, you should find which column is trial index
        log_rt = log_rt * 2; % rt in eeg_log(c15) * 2 matches the actual rt(c6) in behavioral data [Pei, May 2020]
        
        % load the information from BVAinfo (the first column is reaction time, second column is triggers)
        rt_afride = infoBVA(:,1);
        % infoBVA2 = infoBVA(:,1); % for later matching with corrected_log [Pei, June 2020] 
        log_ind2 = log_ind;
        
        % normal: up to k
        if (j ~= [7 8 18 19 20 21 22 23])
        for k = 1:length(log_ind)-1
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial 
                end
            end
        end
        end

        % up to k+4
        if (j==7 & i ~= 24) | (j==8 & i ~= [8 9 13 21 22]) ...
                | (j==18 & i ~= [8 14 15 24]) | (j==19 & i ~= [16 21])...
                | (j==20 & i ~= [5 11 19]) | (j==21 & i ~= [9 13 19 21 22 23])...
                | (j==22 & i ~= [8 24])| (j==23 & i ~= [5 9 10 13 19])
        for k = 1:length(log_ind)-1
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) ...% if RTs don't match
                        & (rt_afride(k+1) == log_rt(k)) % BUT the next row matches
                    rt_afride(k) = []; % delete current row
                    log_ind2 = [log_ind2(1:k,:); nan(1); log_ind2(k+1:end,:)]; % insert 1 NaN [Pei, June 2020]
                    
                elseif (rt_afride(k) ~= log_rt(k)) ...
                        & (rt_afride(k+2) == log_rt(k))
                    rt_afride(k:k+1) = [];
                    log_ind2 = [log_ind2(1:k,:); nan(2,1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) ...
                        & (rt_afride(k+3) == log_rt(k))
                    rt_afride(k:k+2) = [];
                    log_ind2 = [log_ind2(1:k,:); nan(3,1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) ...
                        & (rt_afride(k+4) == log_rt(k))
                    rt_afride(k:k+3) = [];
                    log_ind2 = [log_ind2(1:k,:); nan(4,1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial
                end
            end
        end
        end
        
        % up to k+3
        if (j==7 & i == 24) | (j==18 & i == [14 24]) | (j==19 & i == 16)...
                | (j==23 & i == [5 19])
        for k = 1:length(log_ind)-1
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) ...% if RTs don't match
                        & (rt_afride(k+1) == log_rt(k)) % BUT the next row matches
                    rt_afride(k) = []; % delete current row
                    log_ind2 = [log_ind2(1:k,:); nan(1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) ...
                        & (rt_afride(k+2) == log_rt(k))
                    rt_afride(k:k+1) = [];
                    log_ind2 = [log_ind2(1:k,:); nan(2,1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) ...
                        & (rt_afride(k+3) == log_rt(k))
                    rt_afride(k:k+2) = [];
                    log_ind2 = [log_ind2(1:k,:); nan(3,1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial
                end
            end
        end
        end
       
        % up to k+2
        if (j==8 & i == [13 21 22]) | (j==18 & i == 15)...
                | (j==21 & i == [9 13 19 21 22])| (j==23 & i == 13)
        for k = 1:length(log_ind)-1
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) ...% if RTs don't match
                        & (rt_afride(k+1) == log_rt(k)) % BUT the next row matches
                    rt_afride(k) = []; % delete current row
                    log_ind2 = [log_ind2(1:k,:); nan(1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) ...
                        & (rt_afride(k+2) == log_rt(k))
                    rt_afride(k:k+1) = [];
                    log_ind2 = [log_ind2(1:k,:); nan(2,1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial
                end
            end
        end
        end
        
        % up to k+1
        if (j==8 & i == [8 9]) | (j==18 & i == 8)| (j==19 & i == 21)...
                | (j==20 & i == 5)| (j==22 & i == 24)| (j==23 & i == [9 10])
        for k = 1:length(log_ind)-1
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) ...% if RTs don't match
                        & (rt_afride(k+1) == log_rt(k)) % BUT the next row matches
                    rt_afride(k) = []; % delete current row
                    log_ind2 = [log_ind2(1:k,:); nan(1); log_ind2(k+1:end,:)]; 
                    
                elseif (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial
                end
            end
        end
        end
        
        % up to k
        if (j==20 & i == [11 19])| (j==21 & i == 23)| (j==22 & i == 8)
        for k = 1:length(log_ind)-1
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial
                end
            end
        end
        end
        
        for k = length(log_ind) % last row
            if k <= length(rt_afride) & ~isnan(log_rt(k)) % ignore existing NANs [Pei, May 2020]
                if (rt_afride(k) ~= log_rt(k)) % if RTs don't match
                    rt_afride(k) = nan(1);
                    log_rt (k) = nan(1);
                    log_ind(k) = nan(1); % index for each trial   
                end
            end
        end
        
        if size(rt_afride,1) < length(log_ind)
            last_trial = size(rt_afride,1);
            log_ind((last_trial+1) : length(log_ind)) = nan(1);
        end
        
        % correct for deleted rows in rt_afride [Pei, June 2020]
        for k = 1:length(log_ind)
            if ~isnan(log_ind(k))
                if (log_ind(k) ~= log_ind2(k))% if values don't match b/w log_ind & log_ind2
                    log_ind((k+1):(end+1)) = log_ind(k:end); % insert new row for NAN
                    log_ind(k) = 99999;% insert new row marking "delete"
                end
            end
        end
        corrected_log = 1:14;
        for k = 1:length(log_ind)
            if log_ind(k) == 99999 % for deleted rows
                corrected_log(k,1:end) = nan(1,14); % mark as NAN
            elseif ~isnan(log_ind(k)) % correct rows
                corrected_log(k,1:end)= log(log_ind(k),:);
            elseif isnan(log_ind(k)) % incorrect rows 
                corrected_log(k,1:end) = 88888; % skip
            end
        end
        % corrected_log = log(log_ind(~isnan(log_ind)),:);
        rowsToDelete = any(corrected_log==88888, 2);
        corrected_log(rowsToDelete,:) = [];
       
        % get the corrected log, the last column is the reaction time same as the rt for RIDE correction
        corrected_log(:,end+1) = infoBVA(:,1);
                
        final_log = [final_log; corrected_log];        
        clear corrected_log infoBVA cond_ind 
        
    end
    save([ridefolder,'final_log\', sub{i}, '_final_log.mat'], 'final_log');

end

clear final_log_dir i j 

%%% Result: Log-data was sorted by condition (just like EEG-data) and both
%%% have now been matched and new log-file (final_log) was created
%%% including infos from both



%%

for i = 1:length(sub)
    
    load([ridefolder, 'log_c\', sub{i}, '_log.mat']);
    final_log = [];
    
    for j = 1:length(con)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This part need to be change according to different experiemnts 
        % find conditions in the log file (cause they are not sorted)
        
        % factor 1
        gng_count = str2double(regexp(con{j}, '\d$', 'match')) + 60;
        if isempty(gng_count)
            gng_count = [61:65];
        end
        
        % factor 2
        if ismember('S', con{j})
            bed = 201;
        elseif ismember('J', con{j})
            bed = 200;
        else
            bed = [200:202];
        end
        
        % factor 3
        dg = str2double(regexp(con{j}, '(?<=B)\d{1}', 'match'));
        if isempty(dg)
            dg = [1:3];
        end
        
        % this code should catch all the trials in one condition (index the coloumn info for each factor in log)
        cond_ind = find(ismember(log(:,9), gng_count) & ismember(log(:,6), bed) & ismember(log(:,7), dg));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % load the file from BVA_info
        load([sub{i},'_',con{j},'_infoBVA.mat']);
        
        % select the information for matching from log file
        log_ind = log(cond_ind, end);  % extract the trial index for each condition, you should find which column is trial index
        log_rt = log(cond_ind,end-2);  % extract the special trigger for each condition, you should find which column is trial index
        log_rt = log_rt * 2; % rt in eeg_log * 2 matches the actual rt in behavioral data [Pei, May 2020]
        
        % load the information from BVAinfo (the first column is reaction time, second column is triggers)
        rt_afride = infoBVA(:,1);
        
        for k = 1:length(log_ind)
            
            if k<=length(rt_afride) & (rt_afride(k) ~= log_rt (k)) % if RTs don't match
                rt_afride((k+1):(end+1)) = rt_afride(k:end); % insert new row for NAN
                rt_afride(k) = nan(1); % mark as NAN in rt_afride
                log_rt (k) = nan(1); % mark as NAN in log
                log_ind(k) = nan(1); % mark the index of the trial as NAN
            end
       
        end
        
        if size(rt_afride,1) < length(log_ind)
            last_trial = size(rt_afride,1);
            log_ind((last_trial+1) : length(log_ind)) = nan(1);
        end
        
        % get the corrected log, the last column is the reaction time same as the rt for RIDE correction
        corrected_log = log(log_ind(~isnan(log_ind)),:);
        corrected_log(:,end+1) = infoBVA(:,1);
                
        final_log = [final_log; corrected_log];        
        clear corrected_log infoBVA cond_ind 
        
    end
    save([ridefolder,'final_log\', sub{i}, '_final_log.mat'], 'final_log');

end

clear final_log_dir i j 

%%% Result: Log-data was sorted by condition (just like EEG-data) and both
%%% have now been matched and new log-file (final_log) was created
%%% including infos from both
%%%%%% Note: the sum of trials in 'final_log' for each participant 
%%%%%% should NOT extend the totoal trial numbers in your experiment. 
%%%%%% If so, troubleshoot! [Pei, July 2020]

%% Concatenate the log file of all subjects

log_all = [];
for i = 1: length(sub)    
    load([ridefolder, 'final_log\', sub{i}, '_final_log.mat']);
    log_all = [log_all; final_log];       
end
save([ridefolder, 'final_log\log_all.mat'], 'log_all');

clear i

%% Export time window average or time series to access in R for analyses  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATH = ridefolder; % edit to your needs
EXPORTPATH = 'single_trial_small\'; % Export folder (segmented data will be saved there)
bl=250; % length of prestimulus interval in ms

% export mean 350 - 500 ms, you can set your own interested time window
mint=200; % export start in ms
maxt=350; % export end in ms
CompName = 'S'; % give a name for the component you interested, R or S

% these lines for test purpose with single subject only
dataset_name = 'S_all.mat';% set the name of dataset you want to calculate ERP, in this project, S_all or R_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the data
BFS = load(strcat(PATH,EXPORTPATH,dataset_name));
BFS = struct2cell(BFS);
BFS = BFS{:};

% number of channels
nchans = size(BFS,1);

% the rest is automatic
EXPERP=nanmean(BFS(:,(mint+bl)*srate/1000:(maxt+bl)*srate/1000,:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [CompName,sprintf('%d%d',mint,maxt) '=EXPERP']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH, CompName, mint, maxt), sprintf('%s%d%d', CompName, mint, maxt), '-v6');



%% Combine the time window we are interested in with log file
% in this step, we need two folders, one is final_log folder, the other is single_trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert time window file, see Manual 6.1. (this file is located in single_trials folder)
single_trial_filename = 'S200350.mat';
% this file is located in final_log folder
sorted_log_filename = 'log_all.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_dir = fullfile(ridefolder, 'export_R');
if exist(R_dir,'dir')~=7
    mkdir(R_dir)
end

% file we want to load
log = load(fullfile(ridefolder, 'final_log', sorted_log_filename)); 
log = struct2cell(log);
log = log{:};

% the file we want to load
tw_eeg = load(fullfile(ridefolder, 'single_trials', single_trial_filename));
filename = fieldnames(tw_eeg);
filename = char(filename);
tw_eeg = struct2cell(tw_eeg);
tw_eeg = tw_eeg{:};

% concatenate log and EEG
data = [log tw_eeg];

save([ridefolder, 'export_R\log_' filename '.mat'], 'data'); % save them in the final_log folder, example, 'log_S250400'
clear R_dir log tw_eeg filename 

%% Export interested channals and convert mat file to txt file
% in this step, we need the folder, export_R and original log file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('chanlocs.mat');    
chanlabels = {chanlocs.labels}; % get electrode lables (ONLY head EEG)
nchans = length(chanlabels);% number of channels
chan = {'Fp1'}; % Put whatever electrode you want. Example, chan = {'Fpz', 'F8', 'Cz'}; if you want all the electrode, then chan = chanlabels;
sub_ind_pos = 1; % this is for the column has subject information

% headline for the txt file, the last column is the mean of all electrodes
headline = [{'trial', 'subject', 'itemID', 'cat_num', 'cat', 'member_nr', 'pic', 'condition', 'head', 'cat_type', 'item_type',...
    'LSA_head', 'type_freq_head', 'lemma_freq_head', 'neighbors', 'semantic_rating','cat_nr_trig', 'cond_trig', 'block_trig',...
    'ordinal_trig', 'rep_trig', 'lag_trig', 'lag', 'cat_slot', 'ordinal_pos', 'list_trig', 'cat_slot_trig', 'RT', 'error'} chan {'mean'}];

% load the full data with interested time window
filename = 'log_S200350';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(fullfile(ridefolder, 'export_R'));
% load the data need to write to txt file
load(fullfile(ridefolder, 'export_R', filename)); 

% Open the mat file we want to export with interested electrode
log_columns = 1: (size(data,2) - nchans); % column index of the log file
chan_columns = size(data,2) - nchans + find(ismember(chanlabels, chan)); % column index of the channels
data = data(:,[log_columns chan_columns]); % get the information of all interested electrode 
data = [data mean(data(:,(size(data,2) - length(chan) + 1):end),2)]; % TG: this is to get the mean of three electrodes and put them in the data 

if length(chan) < nchans % if we want to get a txt file for all electrode, then we do not need to save the mat file, because the original mat file have all the electrodes
    save([strjoin([filename chan], '_') '.mat'], 'data'); % save this as a new data set, example name 'R250400_Fpz_F8_Cz.mat'
end

if length(chan) ~= nchans
    fid = fopen([strjoin([filename chan], '_') '.txt'], 'a'); % open the txt file for write , put in name you want
else
    fid = fopen([filename '.txt'], 'a');
end

for i = 1:(length(headline)-1)
    fprintf(fid, '%15s', headline{i}); % write the headline into txt file
end

fprintf(fid, '%15s\r\n', headline{end});

% clear useless variables
clear i 

% reading text file in log file with coresponding name 
for sj = 1: length(sub)
    ind = str2double(cell2mat(regexp(sub{sj}, '\d*','match'))); % get the index of subject
    fid1 = fopen([num2str(ind) '.txt'], 'r');
    log_line = fgetl(fid1); % headline
    %fclose(fid);
    
    % sort the log file based on the single trial index, which is always the last second column
    sorted_log = sortrows(data(find(data(:,sub_ind_pos) == ind),:), (size(data,2)-length(chan)-2));
    line_num = 0;
    for line = 1:max(sorted_log(:,(size(data,2)-length(chan)-2)))
        
        log_line = fgetl(fid1); % headline
        log_line = regexp(log_line,'\s','split'); % '\t' = tab; '\s' = space [Pei, May 2020]
        log_line = log_line(~cellfun('isempty',log_line));
        if ismember(line, sorted_log(:, (size(data,2)-length(chan)-2)))
            line_num = line_num + 1;
            for i = 1:length(log_line);
                fprintf(fid, '%15s', log_line{i}); % Write formatted data to text. Check 'doc fprintf'. [Pei, May 2020]
                % e.g.'%4.2f' specifies that the first value in each line of output is a floating-point number 
                % with a field width of four digits, including two digits after the decimal point.
            end            
            fprintf(fid, [repmat('%15.4f',1,(length(chan)+1)), '\r\n'], sorted_log(line_num,(size(data,2)-length(chan)):end));
        end
        
    end
    % close the txt file
    fclose(fid1);
    clear sorted_log line_num fid1 log_line 
 
end

fclose(fid);
clear i sj line
