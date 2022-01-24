%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath(path_eeglab)

% select relevant files, here baseline blocks
files=dir([data_path filesep 'MWMRI*clean.set']);

%% loop across trials for baseline blocks
redo=1;
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    SubID=SubID(1:findstr(SubID,'clean')-1);
    if redo==0 && exist([save_path filesep 'SW_' SubID '.mat'])~=0
        continue;
    end
   
    EEG = pop_loadset('filename',files(nF).name,'filepath',files(nF).folder);
    fprintf('... ... duration according to EEG data %g minutes\n',size(EEG.data,2)/EEG.srate/60)
    
    all_Waves=[];
    temp_data=EEG.data;
    temp_data(match_str({EEG.chanlocs.labels},{'ECG'}),:)=[];
    ChanLabels={EEG.chanlocs.labels};
    ChanLabels(match_str({EEG.chanlocs.labels},{'ECG'}))=[];
    temp_data=temp_data-repmat(mean(EEG.data(match_str(ChanLabels,{'TP7','TP8'}),:),1),size(temp_data,1),1);
    temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
    
    [twa_results]=twalldetectnew_TA_v2(temp_data,EEG.srate,0);
    for nE=1:length(ChanLabels)
        all_Waves=[all_Waves ; [repmat([1 nF nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
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
    fprintf('\n')
    save([save_path filesep 'allSW_' SubID],'all_Waves','ChanLabels')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=75/2;
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./EEG.srate);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_' SubID],'slow_Waves','paramSW','ChanLabels')
    
end

