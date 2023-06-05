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
redo=1;
all_averp_data=[];
nFc=0;
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
        
        % select trials
        new_test_res=[];
        maxlag_in=[]; maxlag_out=[];
        for nP=1:length(probe_events_sample)
            this_probe_time=probe_res(nP,2);
            temp_test_res=test_res(test_res(:,8)<this_probe_time & test_res(:,8)>this_probe_time-25,:);
            temp_test_res(:,2)=nP;
            temp_test_res(:,4)=temp_test_res(:,4)-temp_test_res(end,4);
            temp_test_res(:,8)=temp_test_res(:,8)-this_probe_time;
            
            temp_trials=(trial_events_sample(trial_events_epoch==nP)-probe_events_sample(nP))/EEG.srate;
            for nT=1:length(temp_trials)
                [closestvalue,index]=findclosest(temp_test_res(:,8)',temp_trials(nT));
                if abs(temp_trials(nT)-closestvalue)<0.1
                    new_test_res=[new_test_res ; temp_test_res(index,:)];
                    maxlag_in(end+1)=closestvalue-temp_trials(nT);
                else
                    maxlag_out(end+1)=closestvalue-temp_trials(nT);
                end
            end
        end
        
        trial_events_sample_rel=trial_events_sample;
        for nP=1:length(probe_events_sample)
            trial_events_sample_rel(trial_events_epoch==nP)=trial_events_sample_rel(trial_events_epoch==nP)-probe_events_sample(nP);
        end
        %     erp_data2=data;
        %     erp_data2.trial=cell(0,0);
        %      erp_data2.time=cell(0,0);
        %      erp_data2.sampleinfo=[];
        %    for nT=1:length(trial_events_sample_rel)
        %         if min(trial_events_sample_rel(nT)-data.time{1}(1)*data.fsample+(-1*data.fsample:2*data.fsample))>0 && max(trial_events_sample_rel(nT)-data.time{1}(1)*data.fsample+(-1*data.fsample:2*data.fsample))<length(data.time{1})
        %             erp_data2.trial{end+1}=data.trial{trial_events_epoch(nT)}(:,trial_events_sample_rel(nT)-data.time{1}(1)*data.fsample+(-1*data.fsample:2*data.fsample));
        %             erp_data2.time{end+1}=-1:1/data.fsample:2;
        %             erp_data2.sampleinfo(end+1,:)=[min(trial_events_sample_rel(nT)-data.time{1}(1)*data.fsample+(-1*data.fsample:2*data.fsample)) max(trial_events_sample_rel(nT)-data.time{1}(1)*data.fsample+(-1*data.fsample:2*data.fsample))];
        %        end
        %     end
        %
        cfg=[];
        cfg.trl=[repmat([-0.5 1.5],length(trial_events_sample),1).*EEG.srate+trial_events_sample' -1*ones(length(trial_events_sample),1)*EEG.srate trial_events_epoch'];
        %         cfg.trl=[repmat([-1 2],length(probe_events_sample),1).*EEG.srate+probe_events_sample' -1*ones(length(probe_events_sample),1)*EEG.srate (1:length(probe_events_sample))'];
        erp_data=ft_redefinetrial(cfg,data);
        
        cfg=[]; % currently ref to 32 (F1???)
        cfg.demean        = 'yes';
        cfg.baselinewindow        = [-0.1 0];
        erp_data=ft_preprocessing(cfg,erp_data);
        
        cfgerp=[];
        cfgerp.keeptrials='yes';
        averp_data=ft_timelockanalysis(cfgerp,erp_data);
        save([files(nF).folder filesep 'ERP_' files(nF).name(1:end-4)],'erp_data');
    else
        nFc=nFc+1;
        load([files(nF).folder filesep 'ERP_' files(nF).name(1:end-4)])
        
        cfgerp=[];
        cfgerp.keeptrials='yes';
        averp_data=ft_timelockanalysis(cfgerp,erp_data);
    end
    all_averp_data=cat(3,all_averp_data,squeeze(nanmean(averp_data.trial,1)));
    
    temp_ERP=squeeze(averp_data.trial(:,19,:));
    temp_ERP=temp_ERP(nanmax(abs(temp_ERP),[],2)<100,:);
    all_Pz_data(nFc,:)=nanmean(temp_ERP,1);
    n_Pz_data(nFc)=size(temp_ERP,1);
    %     pause;
end

