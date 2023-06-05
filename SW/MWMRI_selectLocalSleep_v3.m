%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath((path_fieldtrip))
ft_defaults
addpath(genpath(exgauss_path))
addpath(genpath(FMINSEARCHBND_path))

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean.set']);

myERP_Elec={'Fz','Cz','Pz','Oz'};
States={'ON','MW','MB','??'};
%%
SW_table=array2table(zeros(0,15),'VariableNames',{'SubID','Block','Probe','Elec','State','Vigilance','Misses','FAs','Hit_RT','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.Elec=categorical(SW_table.Elec);
SW_table.State=categorical(SW_table.State);


%% loop across trials for baseline blocks
redo=1;
mean_SW_ERP_byElec=[];
nFc=0;
all_SW_probes=[];
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    sep=findstr(SubID,'clean.set');
    if ismember(SubID,{'MWMRI223','MWMRI243'})
        continue;
    end
    if isempty(sep)
        SubID=SubID(1:end-9);
    else
        SubID=SubID(1:sep(1)-1);
    end
    if redo==0 && exist([save_path filesep 'SW_' SubID '.mat'])~=0
        continue;
    end
    if exist([save_path filesep 'newallSW_' SubID '.mat'])==0
        continue;
    end
    
    % load behaviour
    file_behav=dir([data_path filesep '..' filesep '..' filesep 'Behav' filesep 'wanderIM_behavres_s' SubID(6:end) '*.mat']);
    if ~isempty(file_behav)
        load([file_behav.folder filesep file_behav.name])
    else
        continue;
    end
    
    % load EEG
    EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
    if size(EEG.data,1)<64
        continue;
    end
    ChanLabels={EEG.chanlocs.labels};
    EEG.data=EEG.data(~ismember(ChanLabels,{'ECG','ECG2'}),:,:);
    ChanLabels=ChanLabels(~ismember(ChanLabels,{'ECG','ECG2'}));
    %     evt = ft_read_event([files(nF).folder filesep files(nF).name]);
    %     data = ft_read_data([files(nF).folder filesep files(nF).name]);
    fprintf('... ... duration according to EEG data %g minutes\n',size(EEG.data,2)/EEG.srate/60)
    %      size(EEG.data)
    nFc=nFc+1;
    
    load([save_path filesep 'newallSW_' SubID])
    
    %%% Epoch by probe
    evt=EEG.event;
    event_type={evt.type};
    event_time=[evt.latency];
    findProbe_idx=match_str({evt.type},'P  1');
    findProbe_times=event_time(findProbe_idx);
    
    probe_EEG=permute(EEG.data(:,EEG.times<0,:),[3 1 2]);
    
    probe_EEG=probe_EEG-repmat(mean(probe_EEG(:,match_str(ChanLabels,{'TP9','TP10'}),:),2),[1,size(probe_EEG,2),1]);
    probe_EEG=probe_EEG-repmat(mean(probe_EEG,3),[1,1,size(probe_EEG,3)]);
    
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 7]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150; %150
    paramSW.max_posampl=Nan; %75 normally
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    %     all_Waves(all_Waves(:,2)==1,:)=[];
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
        
        if length(temp_p2p)<100
            warning(sprintf('not enough waves (<100) to compute a threshold!!! Sub %s Elec %s\n',SubID,chan_labels{nE}))
            ExGauss=[ExGauss ; nan(1,4)];
            continue;
        end
        %         [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_p2p);
        %         bins=0:0.1:paramSW.art_ampl;
        %         p_exgauss=exgauss_pdf(bins,X);
        %         end_gaussian=2*bins(find(p_exgauss==max(p_exgauss)));
        thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        %         if ~isempty(paramSW.fixThr)
        %             thr_Wave(nE)=paramSW.fixThr;
        %         else
        %             thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        %         end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_' SubID],'slow_Waves','paramSW')
    %     slow_Waves(slow_Waves(:,3)==length(ChanLabels),:)=[];
    
    for nBl=1:size(probe_EEG,1)
        slow_Waves_perE=[];
        duration_probe=abs(EEG.times(1))/1000/60;
        Fs=EEG.srate;
        for nE=1:length(ChanLabels)
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl)/duration_probe nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,5))/Fs)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,13)) thr_Wave(nE) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,9)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,11))]];
        end
        
        table_length=size(SW_table,1);
        SW_table.SubID(table_length+(1:length(ChanLabels)))=repmat({SubID},length(ChanLabels),1);
        SW_table.Elec(table_length+(1:length(ChanLabels)))=ChanLabels;
        SW_table.State(table_length+(1:length(ChanLabels)))=repmat(States(probe_res(nBl,17)),length(ChanLabels),1);
        
        SW_table.Probe(table_length+(1:length(ChanLabels)))=repmat(nBl,length(ChanLabels),1);
        SW_table.Vigilance(table_length+(1:length(ChanLabels)))=repmat(probe_res(nBl,19),length(ChanLabels),1);
        SW_table.Block(table_length+(1:length(ChanLabels)))=repmat(probe_res(nBl,5),length(ChanLabels),1);
        
        SW_table.SW_density(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,1);
        SW_table.SW_amplitude(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,2);
        SW_table.SW_frequency(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,3);
        SW_table.SW_downslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,4);
        SW_table.SW_upslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,5);
        SW_table.SW_threshold(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,6);
        
        temp_test_res=test_res(test_res(:,1)==probe_res(nBl,5) & test_res(:,4)>=probe_res(nBl,7)-25 & test_res(:,4)<probe_res(nBl,7),:);
        FA=nanmean(temp_test_res(:,12)==0);
        Misses=nanmean(temp_test_res(:,13)==0);
        Hit_RT=nanmean(temp_test_res(temp_test_res(:,5)~=3,11)-temp_test_res(temp_test_res(:,5)~=3,8));
        SW_table.Misses(table_length+(1:length(ChanLabels)))=repmat(Misses,length(ChanLabels),1);
        SW_table.FAs(table_length+(1:length(ChanLabels)))=repmat(FA,length(ChanLabels),1);
        SW_table.Hit_RT(table_length+(1:length(ChanLabels)))=repmat(Hit_RT,length(ChanLabels),1);
    end
end


%%
figure;
temp_plot=[];
for nP=1:40
    temp_plot(nP)=mean(SW_table.SW_density(SW_table.Elec=='Cz' & SW_table.Probe==nP));
end
plot(1:40,temp_plot);
hold on; format_fig;
for k=[10.5 20.5 30.5]
    line([1 1]*k,ylim,'Color','k','LineStyle','--')
end


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

%%
figure
topo_plot=[];
    for nCh=1:length(layout.label)-2
        topo_plot(nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh})); %squeeze(mean(mean(SW_dens_perProbe,2),1));
    end
    simpleTopoPlot_ft(topo_plot', layout,'labels',[],0,1);
    
%%
figure;
for k=1:4
    subplot(2,2,k)
    topo_plot=[];
    for nCh=1:length(layout.label)-2
        topo_plot(nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Block==k)); %squeeze(mean(mean(SW_dens_perProbe,2),1));
    end
    simpleTopoPlot_ft(topo_plot', layout,'on',[],0,1);
    colorbar;
end
%%
figure;
for k=1:4
    subplot(2,2,k)
    topo_plot=[];
    for nCh=1:length(layout.label)-2
        topo_plot(nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.State==States{k})); %squeeze(mean(mean(SW_dens_perProbe,2),1));
    end
    simpleTopoPlot_ft(topo_plot', layout,'on',[],0,1);
    colorbar;
end

%%
figure;
topo_plot=[];
for nCh=1:length(layout.label)-2
    topo_plot(nCh)=nanmean(SW_table.SW_downslope(SW_table.Elec==layout.label{nCh})); %squeeze(mean(mean(SW_dens_perProbe,2),1));
end
simpleTopoPlot_ft(topo_plot', layout,'on',[],0,1);
colorbar;