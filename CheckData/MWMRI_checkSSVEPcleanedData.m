%%
clear all
close all

run ../localdef.m

%% Add fieldtrip
addpath(path_fieldtrip);
ft_defaults;

%% Import data from EEG lab
% EEG     =    load([data_path filesep 'EEG' filesep 'Preprocessed' filesep 'MWMRI202segfastrbcgfilt1sec.set'], '-mat');
EEG     =    load([data_path filesep 'EEG' filesep 'Preprocessed' filesep 'MWMRI202postICA2.set'], '-mat');

%% Compute basic power spectrum
ch_idx  =    match_str({EEG.chanlocs.labels},{'Oz','O1','O2','POz','PO3','PO4'});
% ref_idx =    match_str({EEG.chanlocs.labels},{'TP9','TP10'});
ref_idx =    1:length({EEG.chanlocs.labels});

ch_EEG=[];
allch_EEG=[];
if ndims(EEG.data)==3
    for k=1:length(ch_idx)
        temp_eeg  =    squeeze(EEG.data(ch_idx(k),:,:)-mean(EEG.data(ref_idx,:,:),1));
        ch_EEG(k,:)  =    reshape(temp_eeg,1,numel(temp_eeg));
    end
    for k=1:length({EEG.chanlocs.labels})
        temp_eeg  =    squeeze(EEG.data(k,:,:)-mean(EEG.data(ref_idx,:,:),1));
        allch_EEG(k,:)  =    reshape(temp_eeg,1,numel(temp_eeg));
    end
else
    ch_EEG  =    squeeze(EEG.data(ch_idx,:)-mean(EEG.data(ref_idx,:),1));
    allch_EEG  =    squeeze(EEG.data(:,:)-mean(EEG.data(ref_idx,:),1));
end

%% Gather events
evt_types   =     {EEG.event.type};
evt_samples =     [EEG.event.latency];

evt_blocks_idx      =    match_str(evt_types,'B  1');
evt_blocks_samples  =    evt_samples(evt_blocks_idx);

evt_probes_idx      =    match_str(evt_types,'P  1');
evt_probes_samples  =    evt_samples(evt_probes_idx);

evt_trials_idx      =    match_str(evt_types,'T  1');
evt_trials_samples  =    evt_samples(evt_trials_idx);

evt_MR_idx          =    match_str(evt_types,'S  1');
evt_MR_samples      =    evt_samples(evt_MR_idx);

evt_probes_idx(evt_probes_idx<evt_blocks_idx(1))        = [];
evt_probes_samples(evt_probes_idx<evt_blocks_idx(1))    = [];
evt_trials_idx(evt_trials_idx<evt_blocks_idx(1))        = [];
evt_trials_samples(evt_trials_idx<evt_blocks_idx(1))    = [];
evt_blocks_samples(evt_blocks_idx<=evt_blocks_idx(1))    = [];
evt_blocks_idx(evt_blocks_idx<=evt_blocks_idx(1))        = [];

fprintf('... Blocks: %g | Probes: %g | Trials: %g\n',length(evt_blocks_samples),length(evt_probes_samples),length(evt_trials_idx))

%% Log SNR and Power
aftProbes_logSNR = [];
befProbes_logSNR = [];
aftProbes_logPow = [];
befProbes_logPow = [];
befProbes_logSNR_allCh = [];
befProbes_logPow_allCh =[ ];
param=[];
param.method='welch';
param.w_window=10*EEG.srate;
param.w_overlap=0.5*param.w_window;
param.w_df=0.01;
param.mindist=0.2;

my_events_samples=round(evt_probes_samples);
fprintf('%2.0f/%2.0f\n',0,length(my_events_samples))
for nProbe=1:length(my_events_samples)
    fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nProbe,length(my_events_samples))
    temp_eeg                =   ch_EEG(1,((my_events_samples(nProbe)-40*EEG.srate):my_events_samples(nProbe)));
    [logSNR, faxis, logpow] =   get_logSNR(temp_eeg,EEG.srate,param);
    befProbes_logSNR        =   [befProbes_logSNR ; logSNR];
    befProbes_logPow        =   [befProbes_logPow ; logpow];
    
    temp_eeg                =   ch_EEG(1,((my_events_samples(nProbe)):(my_events_samples(nProbe)+40*EEG.srate)));
    [logSNR, faxis, logpow] =   get_logSNR(temp_eeg,EEG.srate,param);
    aftProbes_logSNR        =   [aftProbes_logSNR ; logSNR];
    aftProbes_logPow        =   [aftProbes_logPow ; logpow];
    
    for k=1:size(allch_EEG,1)
        temp_eeg                =   allch_EEG(k,((my_events_samples(nProbe)-40*EEG.srate):my_events_samples(nProbe)));
        [logSNR, faxis, logpow] =   get_logSNR(temp_eeg,EEG.srate,param);
        befProbes_logSNR_allCh(nProbe,k,:)  =   logSNR;
        befProbes_logPow_allCh(nProbe,k,:)  =   logpow;
    end
end


figure;
format_fig;
plot(faxis,mean(befProbes_logPow)); xlim([2 30])
% hold on
% plot(faxis,mean(aftProbes_logPow)); xlim([2 30])
xlabel('Freq (Hz)')
ylabel('Log Power')
title('Power')

figure;
format_fig;
plot(faxis,mean(befProbes_logSNR)); xlim([2 30])
% hold on
% plot(faxis,mean(aftProbes_logPow)); xlim([2 30])
xlabel('Freq (Hz)')
ylabel('Log SNR')
title('SNR')

% subplot(1,2,2);
% plot(faxis,mean(befProbes_logSNR)); xlim([2 30])
% hold on
% plot(faxis,mean(aftProbes_logSNR)); xlim([2 30])

%%
[~,F_13_45_idx]=findclosest(faxis,14.86);

cfg = [];
ori_Channels={EEG.chanlocs.labels}; ori_Channels{1}='Fp1';
cfg.layout = 'biosemi64.lay';
cfg.channel=ori_Channels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(mean(mean(befProbes_logPow_allCh(:,match_str(ori_Channels,layout.label{nCh}),faxis>10 & faxis<12),1),3));
end
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
simpleTopoPlot_ft(temp_topo', layout,'on',[],1,1);
colorbar;
% caxis([-1 1]*0.6)
%% ERP
ERP_all = cell(1,size(ch_EEG,1));
ERP_all2 = [];

my_events_samples=round(evt_trials_samples);
     fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',0,length(my_events_samples))
for nTrials=1:length(my_events_samples)
     fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',nTrials,length(my_events_samples))
   for k=1:size(ch_EEG,1)
        temp_eeg       =   ch_EEG(k,((my_events_samples(nTrials)-0.5*EEG.srate):(my_events_samples(nTrials)+1.5*EEG.srate)));
        temp_eeg       =   temp_eeg - mean(temp_eeg(0.3*EEG.srate:0.5*EEG.srate));
        if max(abs(temp_eeg))>80
            continue;
        end
        ERP_all{k}        =   [ERP_all{k} ; temp_eeg];
    end
    
    for k=1:size(allch_EEG,1)
        temp_eeg       =   allch_EEG(k,((my_events_samples(nTrials)-0.5*EEG.srate):(my_events_samples(nTrials)+1.5*EEG.srate)));
        temp_eeg       =   temp_eeg - mean(temp_eeg(0.3*EEG.srate:0.5*EEG.srate));
        if max(abs(temp_eeg))>80
            continue;
        end
        ERP_all2(nTrials,k,:)   =   temp_eeg;
    end
end

figure; hold on;
xTime=-0.5:1/EEG.srate:1.5;
    for k=1:3 %size(ch_EEG,1)
plot(xTime,mean(ERP_all{k}));
    end
   legend({'Oz','O1','O2'}) 
    format_fig;
    xlabel('Time (s)')
    ylabel('ERP (\muV)')

    
%     figure; hold on;
%     for k=6 %size(ch_EEG,1)
%         plot(-0.1:1/EEG.srate:0.5,mean(ERP_all2{k}));
%     end
%     xlim([-0.05 0.3])
    

temp_topo=[];
for nCh=1:length(layout.label)-2
    temp_topo(nCh)=squeeze(mean(mean(ERP_all2(:,match_str(ori_Channels,layout.label{nCh}),xTime>0.4 & xTime<0.6),1),3));
end
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
simpleTopoPlot_ft(temp_topo', layout,'on',[],1,1);
        caxis([-1 1]*1.5)
%%
figure
plot(ch_EEG(1,:));
hold on;
scatter(round(evt_blocks_samples),ch_EEG(1,round(evt_blocks_samples)),'og');
scatter(round(evt_probes_samples),ch_EEG(1,round(evt_probes_samples)),'or');
scatter(round(evt_trials_samples),ch_EEG(1,round(evt_trials_samples)),'ok');
