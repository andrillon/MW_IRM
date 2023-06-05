%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath(genpath(path_eeglab))
addpath((path_fieldtrip))
ft_defaults;
% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean.set']);
behav_path=[data_path filesep '..' filesep '..' filesep 'Behav'];
behavfiles=dir([behav_path filesep filesep 'wanderIM_behavres_s*.mat']);

myERP_Elec={'Fz','Cz','Pz','Oz'};

%% loop across trials for baseline blocks
redo=0;
all_averp_data=[];
nFc=0;
for nF=3:length(files)
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
    
    if redo==1 || exist([files(nF).folder filesep 'ERP_' files(nF).name(1:end-4) '.mat'])==0
        nFc=nFc+1;
        EEG = pop_loadset('filename',files(nF).name,'filepath',files(nF).folder);
        fprintf('... ... dimension of EEG data %g - %g - %g\n',size(EEG.data,1),size(EEG.data,2),size(EEG.data,3))
        fprintf('... ... duration according to EEG data %g minutes\n',size(EEG.data,2)/EEG.srate/60)
        data = eeglab2fieldtrip(EEG,'preprocessing','none');
        
        %%% load behav files
        behavfiles_names={behavfiles.name};
        this_file=find_trials(behavfiles_names,SubID(6:8));
        if isempty(this_file)
            warning('PROBLEM BEHAV FILES')
            continue
        end
        load([behav_path filesep behavfiles_names{this_file}])
        %     % epoch alongs probes
        %     cfg=[];
        %     cfg.trl=[];
        %     evt=EEG.event;
        %     event_type={evt.type};
        %     event_time=[evt.latency];
        %     findProbe_idx=match_str({evt.type},'P  1');
        %     findProbe_times=event_time(findProbe_idx);
        %     for nP=1:length(findProbe_idx)
        %         if min((-20*EEG.srate:0*EEG.srate)+round(findProbe_times(nP)))>1
        %             cfg.trl=[cfg.trl ; [-20*EEG.srate+findProbe_times(nP) findProbe_times(nP) -20*EEG.srate]];
        %         end
        %     end
        
        ChanLabels={EEG.chanlocs.labels};
        cfg=[]; % currently ref to 32 (F1???)
        cfg.reref         = 'yes';
        cfg.refchannel ={'TP10','TP10'};
        cfg.channel =ChanLabels(~ismember(ChanLabels,{'ECG','ECG2'}));
        cfg.demean        = 'yes';
        cfg.lpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord      = 4;
        cfg.lpfreq         = 40;
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.hpfilttype     = 'but';
        cfg.hpfiltord      = 4;
        cfg.hpfreq         = 0.1;
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50 100 150];
        data=ft_preprocessing(cfg,data);
        maxAbs=[];
        for nTr=1:length(data.trial)
            maxAbs(nTr,:)=max(abs(data.trial{nTr}),[],2);
        end
        
        % re-define trials
        trial_events_idx=find(ismember({EEG.event.type},'T  1'));
        probe_events_idx=find(ismember({EEG.event.type},'P  1'));
        evt_lat=[EEG.event.latency];
        evt_epoch=[EEG.event.epoch];
        trial_events_sample=evt_lat(trial_events_idx);
        probe_events_sample=evt_lat(probe_events_idx);
        trial_events_epoch=evt_epoch(trial_events_idx);
        %         if length(probe_events_sample)~=size(probe_res,1)
        %             if  length(probe_events_sample)==size(probe_res,1)-1
        %                 probe_res(1,:)=[];
        %             else
        %             warning('PROBLEM NUMBER OF PROBES')
        %             continue;
        %             end
        %         end
        %
        % select trials
        erp_info=[]; out_info=[]; flag_keep=[];
        maxlag_in=[]; maxlag_out=[];
        nPc=1;
        nTc=0;
        for nP=1:length(probe_events_sample)
            temp_trial_events_sample=trial_events_sample(trial_events_epoch==nP);
            for nT=1:length(temp_trial_events_sample)
                this_probe_time=probe_res(nPc,2);
                temp_test_res=test_res(test_res(:,8)<this_probe_time & test_res(:,8)>this_probe_time-25,:);
                while isempty(temp_test_res)
                    nPc=nPc+1;
                    this_probe_time=probe_res(nPc,2);
                    temp_test_res=test_res(test_res(:,8)<this_probe_time & test_res(:,8)>this_probe_time-25,:);
                end
                temp_test_res(:,2)=nPc;
                temp_test_res(:,4)=temp_test_res(:,4)-temp_test_res(end,4);
                temp_test_res(:,8)=temp_test_res(:,8)-this_probe_time;
                
                temp_trial=(temp_trial_events_sample(nT)-probe_events_sample(nP))/EEG.srate;
                [closestvalue,index]=findclosest(temp_test_res(:,8)',temp_trial);
                nTc=nTc+1;
                if abs(temp_trial-closestvalue)<0.1
                    erp_info=[erp_info ; [temp_test_res(index,:) probe_res(nPc,:)]];
                    maxlag_in(end+1)=closestvalue-temp_trial;
                    flag_keep(nTc)=1;
                else
                    out_info=[out_info ; [temp_test_res(index,:) probe_res(nPc,:)]];
                    maxlag_out(end+1)=closestvalue-temp_trial;
                    flag_keep(nTc)=0;
                end
            end
            nPc=nPc+1;
        end
        fprintf('... ... %g trials discarded\n',sum(flag_keep==0))
        
        
        cfg=[];
        trial_events_sample=trial_events_sample(flag_keep==1);
        trial_events_epoch=trial_events_epoch(flag_keep==1);
        cfg.trl=[repmat([-0.5 1.5],length(trial_events_sample),1).*EEG.srate+trial_events_sample' -0.5*ones(length(trial_events_sample),1)*EEG.srate trial_events_epoch'];
        %         cfg.trl=[repmat([-1 2],length(probe_events_sample),1).*EEG.srate+probe_events_sample' -1*ones(length(probe_events_sample),1)*EEG.srate (1:length(probe_events_sample))'];
        erp_data=ft_redefinetrial(cfg,data);
        
        cfg=[]; % currently ref to 32 (F1???)
        cfg.demean        = 'yes';
        cfg.baselinewindow        = [-0.1 0];
        erp_data=ft_preprocessing(cfg,erp_data);
        
        
        cfgerp=[];
        cfgerp.keeptrials='yes';
        averp_data=ft_timelockanalysis(cfgerp,erp_data);
        
        if size(averp_data.trial,1)~=size(erp_info,1)
            warning('PROBLEM NUMBER OF TRIALS in ERP/INFO')
            continue;
        end
        
        save([files(nF).folder filesep 'ERP_' files(nF).name(1:end-4)],'averp_data','erp_info');
    else
        nFc=nFc+1;
        load([files(nF).folder filesep 'ERP_' files(nF).name(1:end-4)])
    end
    all_averp_data=cat(3,all_averp_data,squeeze(nanmean(averp_data.trial,1)));
    
    temp_ERP=squeeze(averp_data.trial(:,18,:));
    temp_ERP(nanmax(abs(temp_ERP),[],2)>100,:)=nan;
    all_Pz_data(nFc,:)=nanmean(temp_ERP,1);
    all_Pz_data_byC(nFc,1,:)=nanmean(temp_ERP(erp_info(:,5)==3,:),1);
    all_Pz_data_byC(nFc,2,:)=nanmean(temp_ERP(erp_info(:,5)~=3,:),1);
    n_Pz_data(nFc)=size(temp_ERP,1);
    n_Pz_data_byC(nFc,1)=sum(~isnan(temp_ERP(erp_info(:,5)==3,1)));
    n_Pz_data_byC(nFc,2)=sum(~isnan(temp_ERP(erp_info(:,5)~=3,1)));
    
    all_Pz_data_byC_byS(nFc,1,1,:)=nanmean(temp_ERP(erp_info(:,5)==3 & erp_info(:,30)==1,:),1);
    all_Pz_data_byC_byS(nFc,1,2,:)=nanmean(temp_ERP(erp_info(:,5)==3 & erp_info(:,30)==2,:),1);
    all_Pz_data_byC_byS(nFc,1,3,:)=nanmean(temp_ERP(erp_info(:,5)==3 & erp_info(:,30)==3,:),1);
    
    all_Pz_data_byC_byS(nFc,2,1,:)=nanmean(temp_ERP(erp_info(:,5)~=3 & erp_info(:,30)==1,:),1);
    all_Pz_data_byC_byS(nFc,2,2,:)=nanmean(temp_ERP(erp_info(:,5)~=3 & erp_info(:,30)==2,:),1);
    all_Pz_data_byC_byS(nFc,2,3,:)=nanmean(temp_ERP(erp_info(:,5)~=3 & erp_info(:,30)==3,:),1);
    %     pause;
end

%%
figure
plot(averp_data.time,squeeze(nanmean(all_Pz_data_byC,1))');
legend({'Go','NoGo'});
format_fig;

%%
figure
for nS=1:3
    subplot(1,3,nS)
plot(averp_data.time,squeeze(nanmean(all_Pz_data_byC_byS(:,:,nS,:),1))');
legend({'Go','NoGo'});
format_fig;
end