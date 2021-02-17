%% clear all
close all
clear all

%% Paths
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
path_data='/Volumes/GoogleDrive/Shared drives/MW_fMRI_EEG/Data/';

addpath(path_fieldtrip);
ft_defaults;

%% Load EEG and behavioural data
name_file='MWMRI202';
sub_ID='202';

hdr=ft_read_header([path_data filesep 'EEG' filesep name_file '.eeg']);
evt=ft_read_event([path_data filesep 'EEG' filesep name_file '.eeg']);

behav_files=dir([path_data filesep 'Behav' filesep 'wanderIM_behavres_s' sub_ID '*.mat']);
load([behav_files.folder filesep behav_files.name])

if name_file=='MWMRI202'
    allBlocks=find(ismember({evt.type},{'B'}));
    evt(1:allBlocks(2)-1)=[];
    probe_res(1:2,:)=[];
    test_res(1:47,:)=[];
end

%% Count the number of events
% TA: Not sure if correct
% A: begining of task
% E: end of task
% B: start of block (expected 4)
% C: end of block (expected 4)
% P: start of a probe
% Q: end of a probe
% S: MRI triggers
% T: start of a trial

fprintf('... found %g start blocks (expected %g)\n',sum(ismember({evt.type},{'B'})),length(unique(test_res(:,1))))
fprintf('... found %g end blocks (expected %g)\n',sum(ismember({evt.type},{'C'})),length(unique(test_res(:,1))))
fprintf('... found %g start probes (expected %g)\n',sum(ismember({evt.type},{'P'})),size(probe_res,1))
fprintf('... found %g end blocks (expected %g)\n',sum(ismember({evt.type},{'Q'})),size(probe_res,1))
fprintf('... found %g start trials (expected %g)\n',sum(ismember({evt.type},{'T'})),size(test_res,1))
fprintf('... found %g MRI triggers\n',sum(ismember({evt.type},{'S'})))

%% Plot the EEG and PTB events
figure;
hold on;
evt_times=[evt.sample]/hdr.Fs;

evt_B_time=evt_times(ismember({evt.type},{'B'}));
evt_P_time=evt_times(ismember({evt.type},{'P'}));
evt_T_time=evt_times(ismember({evt.type},{'T'}));
evt_S_time=evt_times(ismember({evt.type},{'S'}));
evt_times=evt_times-evt_B_time(1);
evt_P_time=evt_P_time-evt_B_time(1);
evt_T_time=evt_T_time-evt_B_time(1);
evt_B_time=evt_B_time-evt_B_time(1);

behav_B_time=all_tstartblock;
behav_P_time=probe_res(:,3);
behav_T_time=test_res(:,8);
behav_B_time=(behav_B_time-all_tstartblock(1));
behav_P_time=(behav_P_time-all_tstartblock(1));
behav_T_time=(behav_T_time-all_tstartblock(1));


plot(evt_times,ismember({evt.type},{'A'}),'k')
plot(evt_times,ismember({evt.type},{'E'}),'k:')

plot(evt_times,ismember({evt.type},{'B'}),'r')
plot(evt_times,ismember({evt.type},{'C'}),'r:')

plot(evt_times,ismember({evt.type},{'P'}),'g')
plot(evt_times,ismember({evt.type},{'Q'}),'g:')

scatter(behav_B_time,ones(1,length(behav_B_time)),'or')
scatter(behav_P_time,ones(1,length(behav_P_time)),'og')
% scatter(behav_T_time,ones(1,length(behav_T_time)),'or')

%% Check distance between PTB and EEG events
diff_Probes=[];
for nPr=1:size(probe_res,1)
    [value, index]=findclosest(evt_P_time,behav_P_time(nPr));
    diff_Probes(nPr)=behav_P_time(nPr)-evt_P_time(index);
end

diff_Trials=[];
for nTr=1:size(test_res,1)
    [value, index]=findclosest(evt_T_time,behav_T_time(nTr));
    diff_Trials(nTr)=behav_T_time(nTr)-evt_T_time(index);
end

figure;
hist(diff_Probes,40);

figure;
hist(diff_Trials,2000);

%% Load EEG
data=ft_read_data([path_data filesep 'EEG' filesep name_file '.eeg'],'header',hdr,'chanindx',17);
%% over block 1
evt_times2=[evt.sample];
evt_B_time2=evt_times2(ismember({evt.type},{'B'}));
evt_C_time2=evt_times2(ismember({evt.type},{'C'}));
figure;
for nBl=1:4
    subplot(2,2,nBl)
    [faxis,pow]=get_PowerSpec(data(1,evt_B_time2(nBl):evt_C_time2(nBl)),hdr.Fs,1,0);
    plot(faxis,pow);
    xlim([1 30]); hold on; for k=1:30, line([1 1]*k*0.6406,ylim,'Color','r','LineStyle','--'); end; xlim([5 20])
end

%%
figure;
for nBl=1:3
    subplot(1,3,nBl)
    [faxis,pow]=get_PowerSpec(data(1,evt_C_time2(nBl):evt_B_time2(nBl+1)),hdr.Fs,1,0);
    plot(faxis,pow);
    xlim([1 30]); hold on; for k=1:30, line([1 1]*k*0.6406,ylim,'Color','r','LineStyle','--'); end; xlim([5 20])
end

figure;
    [faxis,pow]=get_PowerSpec(data(1,evt_C_time2(4):end),hdr.Fs,1,0);
    plot(faxis,pow);
    xlim([1 30]); hold on; for k=1:30, line([1 1]*k*0.6406,ylim,'Color','r','LineStyle','--'); end; xlim([5 20])
