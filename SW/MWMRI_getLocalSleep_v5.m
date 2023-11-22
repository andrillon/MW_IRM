%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% run /Users/kuszti/Documents/Git/MW_IRM/SW/localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath(genpath(path_eeglab))

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean3.set']);

%% loop across trials for baseline blocks
redo=1;
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    sep=findstr(SubID,'clean.set');
    if isempty(sep)
        SubID=SubID(1:end-10); %change hardcoding 
    else
        SubID=SubID(1:sep(1)-1);
    end
%         if ~ismember(SubID,{'MWMRI227'})
%         continue;
%         end
    
    if redo==0 && exist([save_path filesep 'RPAcorr_allSW_' SubID '.mat'])~=0
        continue;
    end
    
    EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
    if isempty(EEG.chanlocs)
        continue;
    end
%     if size(EEG.data,1)<63
%         continue;
%     end
    ChanLabels={EEG.chanlocs.labels};
    EEG.data=EEG.data(~ismember(ChanLabels,{'ECG','ECG2'}),:,:);
    ChanLabels=ChanLabels(~ismember(ChanLabels,{'ECG','ECG2'}));
    %     evt = ft_read_event([files(nF).folder filesep files(nF).name]);
    %     data = ft_read_data([files(nF).folder filesep files(nF).name]);
    fprintf('... ... duration according to EEG data %g minutes - %g EEG channels\n',size(EEG.data,2)/EEG.srate/60,size(EEG.data,1))
%      size(EEG.data)
     
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
    
%     filt_probe_EEG=probe_EEG;
%     for nP=1:40
%         for nC=1:62
%             filt_probe_EEG(nP,nC,:)=bandpass(double(squeeze(probe_EEG(nP,nC,:))),EEG.srate,0.5,10,2);
%         end
%     end
%     for nP=1:40
%         r_eeg=corr(squeeze(filt_probe_EEG(nP,:,:))');
%     end
%     filt_probe_EEG2=probe_EEG2;
%     for nP=1:60
%         for nC=1:63
%             filt_probe_EEG2(nP,nC,:)=bandpass(double(squeeze(probe_EEG2(nP,nC,:))),EEG2.srate,0.5,10,2);
%         end
%     end
%     for nP=1:60
%         r_eeg2=corr(squeeze(filt_probe_EEG2(nP,:,:))');
%     end
%     
%     figure;
%     subplot(2,1,1);
%     imagesc(r_eeg(IA(1:15),IA(1:15)));
%     set(gca,'xTick',1:15,'XTickLabel',ChanLabels(IA(1:15)));
%     set(gca,'yTick',1:15,'yTickLabel',ChanLabels(IA(1:15)));
%     colorbar; caxis([-1 1])
%     title('EEG/fMRI dataset');
%     format_fig;
%     subplot(2,1,2); imagesc(r_eeg(IB(1:15),IB(1:15)));
%     set(gca,'xTick',1:15,'XTickLabel',ChanLabels2(IB(1:15)));
%     set(gca,'yTick',1:15,'yTickLabel',ChanLabels2(IB(1:15)));
%     colorbar; caxis([-1 1])
%  title('EEG NatComm dataset');
%     format_fig;
    
    all_Waves=[];
    for nP=1:size(probe_EEG,1)
        if sum(~isnan(squeeze(probe_EEG(nP,:,1))))==0
            continue;
        end
           fprintf('probe %2.0f/%2.0f\n',nP,size(probe_EEG,1))
%            if mean(max(abs(squeeze(probe_EEG(nP,:,:))),[],2)>250)>0.2
%                continue;
%            end
     [twa_results]=twalldetectnew_TA_v2(squeeze(probe_EEG(nP,:,:)),EEG.srate,0);
        for nE=1:size(probe_EEG,2)
            all_Waves=[all_Waves ; [repmat([1 nP nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
                cell2mat(twa_results.channels(nE).negzx)' ...
                cell2mat(twa_results.channels(nE).poszx)' ...
                cell2mat(twa_results.channels(nE).wvend)' ...
                cell2mat(twa_results.channels(nE).maxnegpk)' ...
                cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
                cell2mat(twa_results.channels(nE).maxpospk)' ...
                cell2mat(twa_results.channels(nE).maxpospkamp)' ...
                cell2mat(twa_results.channels(nE).mxdnslp)' ...
                cell2mat(twa_results.channels(nE).mxupslp)' ...
                cell2mat(twa_results.channels(nE).maxampwn)' ...
                cell2mat(twa_results.channels(nE).minampwn)' ...
                ]];
        end
    end
    fprintf('\n')
    save([save_path filesep 'RPAcorr_allSW_' SubID],'all_Waves','ChanLabels')
    
%     %     %%% clean detection
%     paramSW.prticle_Thr=90; % 80 or 90 or 95
%     paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
%     paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
%     paramSW.fixThr=[];
%     paramSW.art_ampl=150;
%     paramSW.max_posampl=75;
%     paramSW.max_Freq=7;
%     
%     all_Waves=double(all_Waves);
%     all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./EEG.srate);
%     fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
%     fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
%     fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
%     all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
%     
%     thr_Wave=[];
%     slow_Waves=[];
%     for nE=1:64
%         thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
%         temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
%         
%         if ~isempty(paramSW.fixThr)
%             thr_Wave(nE)=paramSW.fixThr;
%         else
%             thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
%         end
%         slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
%     end
    %     save([save_path filesep 'SW_' SubID],'slow_Waves','paramSW','ChanLabels')
    
end

