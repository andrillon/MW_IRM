%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
% close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath(genpath(exgauss_path))
addpath(genpath(FMINSEARCHBND_path))
% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean5.set']);

myERP_Elec={'Fz','Cz','Pz','Oz','C5','C6'};
myERP_Elec2={{'F7','FT9'},{'F8','FT10'}};

load ../Datasetinfo_10-Mar-2023-14-49-21.mat
%% loop across trials for baseline blocks
mean_SW_ERP_byElec=[];
mean_SW_ERP_byElec2=[];
all_ChanLabels=[];
nFc=0;
all_SW_probes=[];
window_before_probes=20; % in seconds

All_slow_Waves_Task = [];
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    sep=findstr(SubID,'clean5.set');
    
    if isempty(sep)
        SubID=SubID(1:end-9);
    else
        SubID=SubID(1:sep(1)-1);
    end
    if exist([save_path filesep 'DSS_allSW_' SubID '.mat'])==0
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
    addpath(genpath(path_eeglab));
    EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
    %EEG = pop_rmbase( EEG, [-25000 0] ,[]);
%     EEG2 = pop_epoch( EEG, {  'R'  }, [-0.5         1.5], 'epochinfo', 'yes');


%     % apply DSS to clean them
%     addpath(genpath('/Users/thandrillon/Work/local/NoiseTools/'))
%     data=EEG2.data;
%     data([32 65],:,:)=[];
%     data=permute(data,[2 1 3]);
%     c0=nt_cov(data);
%     c1=nt_cov(mean(data,3));
%     [todss,pwr0,pwr1]=nt_dss0(c0,c1);
%     NREMOVE=4;
%     z=nt_mmat(data,todss);
%     data_clean=nt_tsr(data,z(:,1:NREMOVE,:)); %
% 
%     % plot var explained
%     figure(1); clf
%     subplot 141;
%     plot(pwr1./pwr0,'Color','k','LineWidth',3);
%     line([1 1]*NREMOVE+0.5,ylim,'Color','r','LineWidth',3)
%     xlabel('Component')
%     % plot results
%     subplot 142;
%     plot(EEG2.times,mean(data,3),'k'); title('data');
%     hold on; plot(EEG2.times,rms(mean(data,3)'),'Color','r','LineWidth',3)
%     subplot 143;
%     plot(EEG2.times,mean(data_clean,3),'k'); title('recovered 1');
%     hold on; plot(EEG2.times,rms(mean(data_clean,3)'),'Color','r','LineWidth',3)
% 
% 
%     data_probe=EEG.data;
%     data_probe([32 65],:,:)=[];
%     data_probe=permute(data_probe,[2 1 3]);
%     z_probe=nt_mmat(data_probe,todss);
%     data_probe_clean=nt_tsr(data_probe,z_probe(:,1:NREMOVE,:)); %
% 
%     EEG3 = EEG;
%     EEG3.data=permute(data_probe_clean,[2 1 3]);
%     EEG3 = pop_epoch( EEG3, {  'R'  }, [-0.5         1.5], 'epochinfo', 'yes');
% %     EEG3 = pop_rmbase( EEG3, [-500 0] ,[]);
% 
%      subplot 144;
%     plot(EEG3.times,mean(EEG3.data,3),'k'); title('recovered 2');
%     hold on; plot(EEG3.times,rms(mean(EEG3.data,3)),'Color','r','LineWidth',3)
    rmpath(genpath(path_eeglab));


%     if size(EEG.data,1)<64
%         continue;
%     end
    ChanLabels={EEG.chanlocs.labels};
    EEG.data=EEG.data(~ismember(ChanLabels,{'ECG','ECG2'}),:,:);
    ChanLabels=ChanLabels(~ismember(ChanLabels,{'ECG','ECG2'}));
    %     evt = ft_read_event([files(nF).folder filesep files(nF).name]);
    %     data = ft_read_data([files(nF).folder filesep files(nF).name]);
    fprintf('... ... duration according to EEG data %g minutes\n',size(EEG.data,2)/EEG.srate/60)
    %      size(EEG.data)
    nFc=nFc+1;
    
    load([save_path filesep 'DSS_allSW_' SubID ])
    
    %%% Epoch by probe
    evt=EEG.event;
    event_type={evt.type};
    event_time=[evt.latency];
    findProbe_idx=match_str({evt.type},'P  1');
    findProbe_times=event_time(findProbe_idx);
    
    probe_EEG=permute(EEG.data,[3 1 2]);
    %     probe_EEG=[];
    %     for nP=1:length(findProbe_times)
    %         if min((-30*EEG.srate:5*EEG.srate)+round(findProbe_times(nP)))>1
    %             probe_EEG(nP,:,:)=EEG.data(:,(-30*EEG.srate:0*EEG.srate-1)+round(findProbe_times(nP)));
    %         else
    %             probe_EEG(nP,:,:)=nan(size(EEG.data,1),30*EEG.srate);
    %         end
    %     end
    
    probe_EEG=probe_EEG-repmat(mean(probe_EEG(:,match_str(ChanLabels,{'TP9','TP10'}),:),2),[1,size(probe_EEG,2),1]);
    probe_EEG=probe_EEG-repmat(mean(probe_EEG,3),[1,1,size(probe_EEG,3)]);
    
    
%     for nP=1:40
%         if sum(all_Waves(:,2)==nP)==0
%             continue;
%         end
%         for nE=1:length(ChanLabels)
%             all_Waves(all_Waves(:,2)==nP & all_Waves(:,3)==nE,16)=mean(squeeze(probe_EEG(nP,:,all_Waves(all_Waves(:,2)==nP & all_Waves(:,3)==nE,8)))>0)';
%         end
%     end
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 7]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150; %150
    paramSW.max_posampl=75; %originally 75 as per the NatCom paper
    paramSW.max_Freq=7;
%     paramSW.min_pptionNeg=1;
    
    all_Waves=double(all_Waves);
    % % Create the histogram
    % figure; % Creates a new figure window
    % histogram(all_Waves(:, 4));
    % xlabel('Peak2Peak Amplitude'); % Label for the x-axis
    % ylabel('Count'); % Label for the y-axis
    % title([SubID, ' Task - Before selection']); % Set the title of the histogram
    % format_fig;
    % % Save the figure
    % outputDir = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/';
    % filename = sprintf('%s%s_SWhist_task_beforeselection.png', outputDir, SubID);  % Construct filename
    % saveas(gcf, filename);  % Save the figure as a PNG file
    % close

    all_Waves(EEG.times(all_Waves(:,5))<-window_before_probes*1000 | EEG.times(all_Waves(:,5))>0,:)=[];
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./EEG.srate);
    fprintf('... ... %g %% waves discarded because of timing\n',mean(all_Waves(:,7)/EEG.srate>30)*100)
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
%     fprintf('... ... %g %% waves discarded because of pption neg elect\n',mean(all_Waves(:,16)>paramSW.min_pptionNeg)*100)
    all_Waves(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
%     all_Waves(all_Waves(:,16)>paramSW.min_pptionNeg | all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
% % Create the histogram
%     figure; % Creates a new figure window
%     histogram(all_Waves(:, 4));
%     xlabel('Peak2Peak Amplitude'); % Label for the x-axis
%     ylabel('Count'); % Label for the y-axis
%     title([SubID, ' Task - After selection']); % Set the title of the histogram
%     format_fig;
%     % Save the figure
%     outputDir = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/';
%     filename = sprintf('%s%s_SWhist_task_afterselection.png', outputDir, SubID);  % Construct filename
%     saveas(gcf, filename);  % Save the figure as a PNG file
%     close
    


% find missing probes
    this_row=find_trials({Dataset.name},SubID);
    if ~isempty(this_row)
        if isempty(Dataset(this_row).BadEpochs)
            probes_missing=[];
        else
            probes_missing= cell2mat(Dataset(this_row).BadEpochs);
        end
    end
    all_probes=1:40;
    all_probes=setdiff(all_probes,probes_missing);

    thr_Wave=[];
    slow_Waves=[];

    for nE=1:length(ChanLabels)
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
     
        thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        %         if ~isempty(paramSW.fixThr)
        %             thr_Wave(nE)=paramSW.fixThr;
        %         else
        %             thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        %         end
   
         slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end

    % % Create the histogram
    % figure; % Creates a new figure window
    % histogram(slow_Waves(:, 4));
    % xlabel('Peak2Peak Amplitude'); % Label for the x-axis
    % ylabel('Count'); % Label for the y-axis
    % title([SubID, ' Task - After thresholding']); % Set the title of the histogram
    % format_fig;
    % % Save the figure
    % outputDir = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/';
    % filename = sprintf('%s%s_SWhist_task_afterthreshold.png', outputDir, SubID);  % Construct filename
    % saveas(gcf, filename);  % Save the figure as a PNG file
    % close

    oldBlockNumbers=slow_Waves(:,2);
    uniqueBlocksSW=unique(oldBlockNumbers);
    for nBl=1:length(uniqueBlocksSW)
        slow_Waves(oldBlockNumbers==uniqueBlocksSW(nBl),2)=all_probes(nBl);
    end
    slow_Waves(:,end+1)=oldBlockNumbers;
    %save([save_path filesep 'prct_25s_DSS_SW_' SubID],'slow_Waves','paramSW','ChanLabels')

    All_slow_Waves_Task = [All_slow_Waves_Task; str2num(SubID(6:end))*ones(size(slow_Waves,1),1), slow_Waves];
    
    [nout,xout]=hist(slow_Waves(:,3),1:length(ChanLabels));
    all_ChanLabels=[all_ChanLabels ; ChanLabels];
    SW_dens(nFc,:)=nout/(window_before_probes/60*40);
    allthr_Wave(nFc,:)=thr_Wave;
    [nout,xout]=hist(all_Waves(:,3),1:length(ChanLabels));
    allSW_dens(nFc,:)=nout/(window_before_probes/60*40);
    
    SW_dens_perProbe(nFc,:,:)=nan(40,size(EEG.data,1));

    SW_dens_perProbe(nFc,1:40,1:63)=nan;
    for nP=all_probes
        [nout,xout]=hist(slow_Waves(slow_Waves(:,2)==nP,3),1:length(ChanLabels));
        SW_dens_perProbe(nFc,nP,:)=nout/(window_before_probes/60);
    end
    
    %%%% check waveform
    
    temp_ERP=cell(2,length(myERP_Elec));
    temp_ERP2=cell(2,length(myERP_Elec2));
    for nE=1:length(ChanLabels)
        all_P2P(nFc,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,4));
        all_NegSl(nFc,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,12));
        all_PosSl(nFc,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,13));
        
        temp_SW=slow_Waves(slow_Waves(:,3)==nE,5);
        temp_SW_nProbe=slow_Waves(slow_Waves(:,3)==nE,16);
        time_SW=EEG.times(temp_SW)/1000;
        
        if ismember(ChanLabels{nE},myERP_Elec)
            % get ERP for slow waves
            for m=1:length(temp_SW)
                if temp_SW(m)-0.5*EEG.srate>0 && temp_SW(m)+1*EEG.srate<25*EEG.srate
                    vec=squeeze(probe_EEG(temp_SW_nProbe(m),nE,(temp_SW(m)-0.5*EEG.srate):(temp_SW(m)+1*EEG.srate)))';
                    vec=vec-mean(vec(1:0.5*EEG.srate));
                    if max(abs(vec))<150
                        temp_ERP{1,find(ismember(myERP_Elec,ChanLabels{nE}))}=[temp_ERP{1,find(ismember(myERP_Elec,ChanLabels{nE}))} ; vec];
                    end
                    
                    vec=squeeze(mean(probe_EEG(temp_SW_nProbe(m),:,(temp_SW(m)-0.5*EEG.srate):(temp_SW(m)+1*EEG.srate)),2))';
                    vec=vec-mean(vec(1:0.5*EEG.srate));
                    if max(abs(vec))<150
                        temp_ERP{2,find(ismember(myERP_Elec,ChanLabels{nE}))}=[temp_ERP{2,find(ismember(myERP_Elec,ChanLabels{nE}))} ; vec];
                    end
                end
            end
        end
        if ismember(ChanLabels{nE},myERP_Elec2{1})
            % get ERP for slow waves
            for m=1:length(temp_SW)
                if temp_SW(m)-0.5*EEG.srate>0 && temp_SW(m)+1*EEG.srate<25*EEG.srate
                    vec=squeeze(probe_EEG(temp_SW_nProbe(m),nE,(temp_SW(m)-0.5*EEG.srate):(temp_SW(m)+1*EEG.srate)))';
                    vec=vec-mean(vec(1:0.5*EEG.srate));
                    if max(abs(vec))<150
                        temp_ERP2{1,find(ismember(myERP_Elec2{1},ChanLabels{nE}))}=[temp_ERP2{1,find(ismember(myERP_Elec2{1},ChanLabels{nE}))} ; vec];
                    end
                end
            end
            this_PairE=match_str(ChanLabels,myERP_Elec2{2}(ismember(myERP_Elec2{1},ChanLabels{nE})));
            temp_SW=slow_Waves(slow_Waves(:,3)==this_PairE,5);
            temp_SW_nProbe=slow_Waves(slow_Waves(:,3)==this_PairE ,16);
            for m=1:length(temp_SW)
                if temp_SW(m)-0.5*EEG.srate>0 && temp_SW(m)+1*EEG.srate<25*EEG.srate
                    vec=squeeze(probe_EEG(temp_SW_nProbe(m),nE,(temp_SW(m)-0.5*EEG.srate):(temp_SW(m)+1*EEG.srate)))';
                    vec=vec-mean(vec(1:0.5*EEG.srate));
                    if max(abs(vec))<150
                        temp_ERP2{2,find(ismember(myERP_Elec2{1},ChanLabels{nE}))}=[temp_ERP2{2,find(ismember(myERP_Elec2{1},ChanLabels{nE}))} ; vec];
                    end
                end
            end
            
            
        end
    end
    temp_SW_ERP_byElec=[];
    for j=1:length(myERP_Elec)
        if ~isempty(temp_ERP{1,j})
            temp_SW_ERP_byElec(j,1,:)=mean(temp_ERP{1,j},1);
            temp_SW_ERP_byElec(j,2,:)=mean(temp_ERP{2,j},1);
        else
            temp_SW_ERP_byElec(j,1,:)=nan(1,751);
            temp_SW_ERP_byElec(j,2,:)=nan(1,751);
        end
    end
    mean_SW_ERP_byElec=cat(4,mean_SW_ERP_byElec,temp_SW_ERP_byElec);
    
    temp_SW_ERP_byElec2=[];
    for j=1:length(myERP_Elec2{1})
        for k=1:2
            if ~isempty(temp_ERP2{k,j})
                temp_SW_ERP_byElec2(j,k,:)=mean(temp_ERP2{k,j},1);
            else
                temp_SW_ERP_byElec2(j,k,:)=nan(1,751);
            end
        end
    end
    mean_SW_ERP_byElec2=cat(4,mean_SW_ERP_byElec2,temp_SW_ERP_byElec2);
    
    probe_res(:,1)=probe_res(:,1)+10*(probe_res(:,5)-1);
    all_SW_probes=[all_SW_probes ; str2num(SubID(6:end))*ones(size(probe_res,1)-sum(ismember(probe_res(:,1),probes_missing)),1) probe_res(~ismember(probe_res(:,1),probes_missing),[1 5 17 18 19]) squeeze(nanmean(SW_dens_perProbe(nFc,~ismember(probe_res(:,1),probes_missing),:),3))' squeeze(SW_dens_perProbe(nFc,~ismember(probe_res(:,1),probes_missing),:))];
end

%%
figure;
set(gcf,'Position',[680   371   560   420]);
for nCh=1:length(myERP_Elec)
    subplot(3,2,nCh);
    temp_plot=squeeze(mean_SW_ERP_byElec(nCh,1,:,:));
    temp_plot2=squeeze(mean_SW_ERP_byElec(nCh,2,:,:));
    
    plot(-0.5:1/EEG.srate:1,temp_plot,'Color','k');
    simpleTplot(-0.5:1/EEG.srate:1,temp_plot',0,[1 0.5 0],0,'-',0.5,1,[],1,4);
    simpleTplot(-0.5:1/EEG.srate:1,temp_plot2',0,[0.5 0.5 0.5],0,'-',0.5,1,[],1,4);
    format_fig;
    xlabel('Time from onset (s)')
    ylabel('Voltage (\muV)')
    title(myERP_Elec{nCh})
    % ylim([-10 3])
end

%%
figure;
simpleTplot(1:40,squeeze(mean(SW_dens_perProbe(:,:,match_str({EEG.chanlocs.labels},'Fz')),3)),0,[1 0.5 0],0,'-',0.5,1,[],1,4);
hold on; format_fig;
for k=[10.5 20.5 30.5]
    line([1 1]*k,ylim,'Color','k','LineStyle','--')
end
xlabel('Probes')
ylabel('SW density')

%%
addpath((path_fieldtrip))
ft_defaults;

ChanLabels={EEG.chanlocs.labels};
ChanLabels(find(ismember(ChanLabels,'FPz')))={'Fpz'};
ChanLabels(find(ismember(ChanLabels,'FP1')))={'Fp1'};

cfg = [];
cfg.layout = 'elec1005.lay';
cfg.center      = 'yes';
cfg.channel = ChanLabels(~ismember(ChanLabels,{'TP9','TP10'}));
layout=ft_prepare_layout(cfg);
% layout.label(match_str(layout.label,{'TP9','TP10'}))=[];
correspCh=[];
for nCh=1:length(layout.label)-2
    correspCh(nCh)=match_str(ChanLabels,layout.label{nCh});
end


%%
figure;
for k=1:size(SW_dens_perProbe,1)
     SubID=files(k).name;
    sep=findstr(SubID,'clean5.set');
    
    if isempty(sep)
        SubID=SubID(1:end-9);
    else
        SubID=SubID(1:sep(1)-1);
    end
    
    subplot(5,5,k);
    topo_plot=squeeze(nanmean(nanmean(SW_dens_perProbe(k,:,:),2),1));
    if sum(isnan(topo_plot))==length(topo_plot)
        continue;
    end
    % topo_plot=squeeze(nanmean(SW_dens,1));
    % topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=0;
    simpleTopoPlot_ft(topo_plot(correspCh), layout,'on',[],0,1);
    colorbar;
    title(SubID)
end
%% Topography

figure;
topo_plot=mean(allthr_Wave,1); %squeeze(mean(mean(SW_dens_perProbe,2),1));
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('Threshold of slow waves detection')

figure;
subplot(2,2,1);
topo_plot=squeeze(nanmean(nanmean(SW_dens_perProbe,2),1));
% topo_plot=squeeze(nanmean(SW_dens,1));
% topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=0;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW density')

% figure;
subplot(2,2,2);
topo_plot=squeeze(nanmean(all_P2P,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW amplitude')

subplot(2,2,3);
topo_plot=squeeze(nanmean(all_PosSl,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW positive slope')

subplot(2,2,4);
topo_plot=squeeze(nanmean(all_NegSl,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW negative slope')

%%
figure;
for nB=1:4
    subplot(1,4,nB);
    topo_plot=squeeze(nanmean(nanmean(SW_dens_perProbe(:,1:10+(nB-1)*10,:),2),1));
    %     topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
    simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
    colorbar;
    %     caxis([6.5 9.5])
end

% % mind-state
% figure;
% States={'ON','MW','MB'};
% for nB=1:3
%     subplot(1,3,nB);
% 
% 
% % Extract field names for all SubIDs in results
% subIDFields = fieldnames(results);
% temp1=[];
% % Loop through each SubID field
% for nsub = 1:length(subIDFields)
%     subIDField = subIDFields{nsub};
%     if isfield(results.(subIDField), 'SART_State_MB')
%         % Extract SART_State_1 data
%         %inds = results.(subIDField).SART_State_ON; %OT
%         %inds = results.(subIDField).SART_State_MW; %OT
%         inds = results.(subIDField).SART_State_MB; %OT
% 
% 
%         temp1(:,nsub) = squeeze(nanmean(SW_dens_perProbe(nsub,inds,:),2));
% 
%     end
% 
% end
% 
%     figure
%     topo_plot=squeeze(nanmean(temp1,2));
%     %     topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
%     simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
%     colorbar;
%     caxis([5 7])
%     title(['Average SW density /min' newline 'Mind-wandering'])
%     %title(['Average SW density /min' newline 'On Task'])
%     %title(['Average SW density /min' newline 'Mind-blanking'])
% 
% 
%     %     caxis([6.5 9.5])
% end
% 
% T = table_SW;
% T.SART_State = categorical(T.SART_State);
% 
% %T.SART_State = mergecats(T.SART_State, {'DK'}, 'MB');
% 
% % Initialize a structure to store results
% results = struct;
% 
% % Get unique combinations of SubID and SART_State
% uniqueSubIDs = unique(T.SubID);
% uniqueSARTStates = unique(T.SART_State);
% 
% % Loop through each SubID
% for i = 1:length(uniqueSubIDs)
%     subID = uniqueSubIDs(i);
% 
%     % Find rows for this SubID
%     subIDRows = T.SubID == subID;
% 
%     % Loop through each SART_State for this SubID
%     for j = 1:length(uniqueSARTStates)
%         sartState = uniqueSARTStates(j);
% 
%         % Find rows for this SART_State and SubID
%         stateRows = T.SART_State == sartState & subIDRows;
% 
%         % Extract unique Probe numbers for this combination
%         uniqueProbes = unique(T.Probe(stateRows));
% 
%         % Store results
%         if ~isempty(uniqueProbes)
%             % Use dynamic field names for the structure
%             subIDField = sprintf('SubID_%d', subID);
%             stateField = sprintf('SART_State_%s', sartState);
% 
%             % Check if field exists for SubID
%             if ~isfield(results, subIDField)
%                 results.(subIDField) = struct;
%             end
% 
%             % Store unique probes under the correct SubID and SART_State
%             results.(subIDField).(stateField) = uniqueProbes;
%             end
%             end
%             end
% 
% 
% 
% 

            %%
figure;
topo_plot=mean(SW_dens(1:end-1,:),1); %squeeze(mean(mean(SW_dens_perProbe,2),1));
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW density')

figure;
topo_plot=squeeze(mean(nanmean(SW_dens_perProbe,2),1));
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW density')


%%
figure;
set(gcf,'Position',[680   371   560   420]);
for nCh=1:size(mean_SW_ERP_byElec2,1)
    subplot(1,2,nCh);
    temp_plot=squeeze(mean_SW_ERP_byElec2(nCh,1,:,:));
    temp_plot2=squeeze(mean_SW_ERP_byElec2(nCh,2,:,:));
    
    plot(-0.5:1/EEG.srate:1,temp_plot,'Color','k');
    simpleTplot(-0.5:1/EEG.srate:1,temp_plot',0,[1 0.5 0],0,'-',0.5,1,[],1,4);
    simpleTplot(-0.5:1/EEG.srate:1,temp_plot2',0,[0.5 0.5 0.5],0,'-',0.5,1,[],1,4);
    format_fig;
    xlabel('Time from onset (s)')
    ylabel('Voltage (\muV)')
    title(myERP_Elec2{1}{nCh})
    % ylim([-10 3])
end