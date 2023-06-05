%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath(genpath(exgauss_path))
addpath(genpath(FMINSEARCHBND_path))

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean.set']);

myERP_Elec={'Fz','Cz','Pz','Oz'};

%% loop across trials for baseline blocks
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
    
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 7]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150; %150
    paramSW.max_posampl=NaN; %originally 75 as per the NatCom paper
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
%     all_Waves(all_Waves(:,2)==1,:)=[];
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./EEG.srate);
    fprintf('... ... %g %% waves discarded because of timing\n',mean(all_Waves(:,7)/EEG.srate>30)*100)
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
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
        [X,fVal,exitFlag,solverOutput] = exgauss_fit(temp_p2p);
        bins=0:0.1:paramSW.art_ampl;
        p_exgauss=exgauss_pdf(bins,X);
        end_gaussian=2*bins(find(p_exgauss==max(p_exgauss)));
%         figure;
%         histogram(zscore(log(temp_p2p)),150,'Normalization','probability');
%         hold on
%         line([1 1]*prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr),ylim,'Color','r')
%         line([1 1]*end_gaussian,ylim,'Color','g')
%         line([1 1]*2,ylim,'Color','g')
%         pause;
%         close(gcf);
        thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        %         if ~isempty(paramSW.fixThr)
        %             thr_Wave(nE)=paramSW.fixThr;
        %         else
        %             thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        %         end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_' SubID],'slow_Waves','paramSW')
    
    [nout,xout]=hist(slow_Waves(slow_Waves(:,7)/EEG.srate<30,3),1:length(ChanLabels));
    SW_dens(nFc,:)=nout/(size(EEG.data,2)*size(EEG.data,1)/EEG.srate/60);
    allthr_Wave(nFc,:)=thr_Wave;
    
    SW_dens_perProbe(nFc,:,:)=nan(40,size(EEG.data,1));
    for nP=1:length(findProbe_times)
        [nout,xout]=hist(slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,7)/EEG.srate<30,3),1:length(ChanLabels));
        SW_dens_perProbe(nFc,nP,:)=nout/(size(EEG.data,2)/EEG.srate/60);
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
                if temp_SW(m)-0.5*EEG.srate>0 && temp_SW(m)+1*EEG.srate<25*EEG.srate
                    vec=squeeze(probe_EEG(temp_SW_nProbe(m),nE,(temp_SW(m)-0.5*EEG.srate):(temp_SW(m)+1*EEG.srate)))';
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
        if ~isempty(temp_ERP{j})
        temp_SW_ERP_byElec(j,:)=mean(temp_ERP{j},1);
        else
            temp_SW_ERP_byElec(j,:)=nan(1,751);
        end
    end
    mean_SW_ERP_byElec=cat(3,mean_SW_ERP_byElec,temp_SW_ERP_byElec);
    
    all_SW_probes=[all_SW_probes ; str2num(SubID(6:end))*ones(size(probe_res,1),1) probe_res(:,[1 5 17 18 19]) squeeze(nanmean(SW_dens_perProbe(nFc,:,:),3))' squeeze(SW_dens_perProbe(nFc,:,:))];
end

%%
figure;
set(gcf,'Position',[680   371   560   420]);
for nCh=1:length(myERP_Elec)
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
xlabel('Probes')
ylabel('SW density')

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
title('Threshold of slow waves detection')

figure;
subplot(2,2,1);
topo_plot=squeeze(mean(nanmean(SW_dens_perProbe,2),1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW density')

subplot(2,2,2);
topo_plot=squeeze(mean(all_P2P,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW amplitude')

subplot(2,2,3);
topo_plot=squeeze(mean(all_PosSl,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW positive slope')

subplot(2,2,4);
topo_plot=squeeze(mean(all_NegSl,1));
topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
colorbar;
title('SW negative slope')

%%
figure;
for nB=1:4
    subplot(1,4,nB);
    topo_plot=squeeze(mean(nanmean(SW_dens_perProbe(:,1:10+(nB-1)*10,:),2),1));
%     topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
    simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
    colorbar;
%     caxis([6.5 9.5])
end