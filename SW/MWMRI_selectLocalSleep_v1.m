%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath(genpath(path_eeglab))

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean.set']);

myERP_Elec={'Fz','Cz','Pz','Oz'};

%% loop across trials for baseline blocks
redo=1;
mean_SW_ERP_byElec=[];
nFc=0;
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    sep=findstr(SubID,'_');
    if isempty(sep)
        SubID=SubID(1:end-9);
    else
        SubID=SubID(1:sep(1)-1);
    end
    if redo==0 && exist([save_path filesep 'SW_' SubID '.mat'])~=0
        continue;
    end
    if exist([save_path filesep 'allSW_' SubID '.mat'])==0
        continue;
    end
    nFc=nFc+1;
    EEG = pop_loadset('filename',files(nF).name,'filepath',files(nF).folder);
    fprintf('... ... duration according to EEG data %g minutes\n',size(EEG.data,2)/EEG.srate/60)
    load([save_path filesep 'allSW_' SubID])
    
    %     %%% load behavioural results and keep only waves 20s before probes
    %     behav_files=dir([root_path filesep 'Behav' filesep 'wanderIM_behavres_s' SubID(end-2:end) '*.mat']);
    %     load([behav_files(end).folder filesep behav_files(end).name])
    
    %%% epoch probe
    temp_data=EEG.data(ismember({EEG.chanlocs.labels},ChanLabels),:);
    temp_data=temp_data-repmat(mean(temp_data(match_str(ChanLabels,{'TP9','TP10'}),:),1),[size(temp_data,1),1]);
    event_type={EEG.event.type};
    event_time=[EEG.event.latency];
    findProbe_idx=match_str({EEG.event.type},'P  1');
    findProbe_times=event_time(findProbe_idx);
    probe_EEG=[];
    for nP=1:length(findProbe_times)
        if min((-30*EEG.srate:5*EEG.srate)+round(findProbe_times(nP)))>1
            probe_EEG(nP,:,:)=temp_data(:,(-30*EEG.srate:5*EEG.srate)+round(findProbe_times(nP)));
        else
            probe_EEG(nP,:,:)=nan(size(temp_data,1),35*EEG.srate+1);
        end
    end
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 7]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./EEG.srate);
    fprintf('... ... %g %% waves discarded because of timing\n',mean(all_Waves(:,7)/EEG.srate>30)*100)
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_Waves(:,7)/EEG.srate>30 | all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:length(ChanLabels)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_' SubID],'slow_Waves','paramSW')
    %     slow_Waves(slow_Waves(:,3)==length(ChanLabels),:)=[];
%     figure('Position',[1 377 1394 420]); histogram(slow_Waves(:,3),1:length(ChanLabels)); set(gca,'xtick',1:length(ChanLabels),'xticklabel',{EEG.chanlocs.labels})
    [nout,xout]=hist(slow_Waves(slow_Waves(:,7)/EEG.srate<30,3),1:length(ChanLabels));
    SW_dens(nFc,:)=nout/(size(EEG.data,2)/EEG.srate/60);
    allthr_Wave(nFc,:)=thr_Wave;
    
    for nP=1:length(findProbe_times)
        [nout,xout]=hist(slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,7)/EEG.srate<30,3),1:length(ChanLabels));
        SW_dens_perProbe(nFc,nP,:)=nout/2;
    end
    
    %%%% check waveform
    labels={EEG.chanlocs.labels};
    
    temp_ERP=cell(1,length(myERP_Elec));
    for nE=1:length(ChanLabels)
        all_P2P(nFc,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,5)/EEG.srate<30,4));
        all_NegSl(nFc,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,5)/EEG.srate<30,12));
        all_PosSl(nFc,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,5)/EEG.srate<30,13));
        
        temp_SW=slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,5)/EEG.srate<30,5);
        temp_SW_nProbe=slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,5)/EEG.srate<30,2);
        time_SW=(temp_SW/EEG.srate)/60;
        
        if ismember(labels(nE),myERP_Elec)
            % get ERP for slow waves
            for m=1:length(temp_SW)
                if temp_SW(m)-0.5*EEG.srate>0 && temp_SW(m)+1*EEG.srate<30*EEG.srate
                    vec=probe_EEG(temp_SW_nProbe(m),nE,(temp_SW(m)-0.5*EEG.srate):(temp_SW(m)+1*EEG.srate));
                    vec=vec-mean(vec(1:0.5*EEG.srate));
                    if max(abs(vec))<150
                        temp_ERP{find(ismember(myERP_Elec,labels(nE)))}=[temp_ERP{find(ismember(myERP_Elec,labels(nE)))} ; vec];
                    end
                end
            end
        end
    end
    temp_SW_ERP_byElec=[];
    for j=1:length(myERP_Elec)
        temp_SW_ERP_byElec(j,:)=mean(temp_ERP{j});
    end
    mean_SW_ERP_byElec=cat(3,mean_SW_ERP_byElec,temp_SW_ERP_byElec);
end

%%
figure;
set(gcf,'Position',[680   371   560   420]);
for nCh=1:4
    subplot(2,2,nCh);
    temp_plot=squeeze(mean_SW_ERP_byElec(nCh,:,:));
    
     plot(-0.5:1/EEG.srate:1,temp_plot,'Color','k');
   simpleTplot(-0.5:1/EEG.srate:1,temp_plot',0,[1 0.5 0],0,'-',0.5,1,[],1,4);
    format_fig;
    xlabel('Time from onset (s)')
    ylabel('Voltage (\muV)')
    title(myERP_Elec{nCh})
% ylim([-10 3])
end

%%
figure;
simpleTplot(1:40,squeeze(mean(SW_dens_perProbe(:,:,match_str({EEG.chanlocs.labels},'Cz')),3)),0,[1 0.5 0],0,'-',0.5,1,[],1,4);
hold on; format_fig;
for k=[10.5 20.5 30.5]
line([1 1]*k,ylim,'Color','k','LineStyle','--')
end

%%
addpath(rmpath(path_eeglab))
addpath((path_fieldtrip))
ft_defaults;

%% Topography


% ChanLabels={EEG.chanlocs.labels};
ChanLabels(find(ismember(ChanLabels,'FPz')))={'Fpz'};
ChanLabels(find(ismember(ChanLabels,'FP1')))={'Fp1'};

cfg = [];
cfg.layout = 'elec1005.lay';
cfg.center      = 'yes';
cfg.channel = ChanLabels;
layout=ft_prepare_layout(cfg);
for nCh=1:length(layout.label)-2
correspCh(nCh)=match_str(ChanLabels,layout.label{nCh});
end

figure;
topo_plot=mean(allthr_Wave,1); %squeeze(mean(mean(SW_dens_perProbe,2),1));
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;

figure;
subplot(2,2,1);
topo_plot=squeeze(mean(SW_dens,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;

subplot(2,2,2);
topo_plot=squeeze(mean(all_P2P,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;

subplot(2,2,3);
topo_plot=squeeze(mean(all_PosSl,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;

subplot(2,2,4);
topo_plot=squeeze(mean(all_NegSl,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
