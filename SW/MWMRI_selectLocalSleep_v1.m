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
files=dir([data_path filesep 'MWMRI*clean_a.set']);

myERP_Elec={'Fz','Cz','Pz','Oz'};

%% loop across trials for baseline blocks
redo=1;
mean_SW_ERP_byElec=[];
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    sep=findstr(SubID,'_');
    SubID=SubID(1:sep(1)-1);
    if redo==0 && exist([save_path filesep 'SW_' SubID '.mat'])~=0
        continue;
    end
    
    EEG = pop_loadset('filename',files(nF).name,'filepath',files(nF).folder);
    fprintf('... ... duration according to EEG data %g minutes\n',size(EEG.data,2)/EEG.srate/60)
    load([save_path filesep 'allSW_' SubID])
    
    %%% load behavioural results and keep only waves 20s before probes
    behav_files=dir([root_path filesep 'Behav' filesep 'wanderIM_behavres_s' SubID(end-2:end) '*.mat']);
    load([behav_files(end).folder filesep behav_files(end).name])

    %%% epoch probe
    temp_data=EEG.data;
    temp_data(match_str({EEG.chanlocs.labels},{'ECG'}),:)=[];
    ChanLabels={EEG.chanlocs.labels};
    ChanLabels(match_str({EEG.chanlocs.labels},{'ECG'}))=[];
     event_type={EEG.event.type};
    event_time=[EEG.event.latency];
    findProbe_idx=match_str({EEG.event.type},'P  1');
    findProbe_times=event_time(findProbe_idx);
    probe_EEG=[];
    nPc=0;
    for nP=1:length(findProbe_times)
        if min((-30*EEG.srate:5*EEG.srate)+round(findProbe_times(nP)))>1
            nPc=nPc+1;
        probe_EEG(nPc,:,:)=temp_data(:,(-30*EEG.srate:5*EEG.srate)+round(findProbe_times(nP)));
        end
    end
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=5;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./EEG.srate);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:63
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
    %     slow_Waves(slow_Waves(:,3)==64,:)=[];
    figure('Position',[1 377 1394 420]); histogram(slow_Waves(:,3),1:64); set(gca,'xtick',1:64,'xticklabel',{EEG.chanlocs.labels})
    [nout,xout]=hist(slow_Waves(:,3),1:63);
    SW_dens(nF,:)=nout/(size(EEG.data,2)/EEG.srate/60);
    allthr_Wave(nF,:)=thr_Wave;
    
    %%%% check waveform
    labels={EEG.chanlocs.labels};
    
    temp_ERP=cell(1,length(myERP_Elec));
    for nE=1:63
        temp_P2P(nF,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,4));
        temp_NegSl(nF,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,12));
        temp_PosSl(nF,nE)=nanmean(slow_Waves(slow_Waves(:,3)==nE,13));
        
        temp_SW=slow_Waves(slow_Waves(:,3)==nE,5);
        temp_SW_nProbe=slow_Waves(slow_Waves(:,3)==nE,2);
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
plot(-0.5:1/EEG.srate:1,squeeze((mean_SW_ERP_byElec))');

