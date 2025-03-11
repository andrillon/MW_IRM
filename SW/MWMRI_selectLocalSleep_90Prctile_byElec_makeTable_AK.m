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
addpath(path_fieldtrip)

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean5.set']);


%%
%table_SW=array2table(zeros(0,25),'VariableNames',{'SubID','Block','Probe','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SART_Miss','SART_FA','SART_HitRT','SART_ON','SART_State','SART_Vig','SART_MWcat','SART_StateError','ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROIfrontR', 'ROIfrontMP', 'ROIbackR', 'ROIbackMP'});
table_SW=array2table(zeros(0,19),'VariableNames',{'SubID','Block','Probe','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SART_Miss','SART_FA','SART_HitRT','SART_ON','SART_State','SART_Vig','SART_MWcat','SART_StateError', 'ROIfront', 'ROIback',});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.Elec=categorical(table_SW.Elec);
table_SW.SART_State=categorical(table_SW.SART_State);
table_SW.SART_StateError=categorical(table_SW.SART_StateError);


States={'ON','MW','MB','DK'};
%% loop across trials for baseline blocks
nFc=0;
all_SW_probes=[];
window_before_probes=20; % in seconds
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
    
    nFc=nFc+1;
    % if nFc==1
    %     addpath(genpath(path_eeglab));
    % EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
    % rmpath(genpath(path_eeglab));
    % end
    
%     load([save_path filesep 'DSS_allSW_' SubID ])
    Fs=500;
    
    %%% clean detection
%     paramSW.prticle_Thr=90; % 80 or 90 or 95
%     paramSW.LimFrqW=[1 7]; % [1 4] or [4 10]
%     paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
%     paramSW.fixThr=[];
%     paramSW.art_ampl=150; %150
%     paramSW.max_posampl=75; %originally 75 as per the NatCom paper
%     paramSW.max_Freq=7;
%     %     paramSW.min_pptionNeg=1;
%     
%     all_Waves=double(all_Waves);
%     all_Waves(EEG.times(all_Waves(:,5))<-window_before_probes*1000 | EEG.times(all_Waves(:,5))>0,:)=[];
%     all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs);
%     fprintf('... ... %g %% waves discarded because of timing\n',mean(all_Waves(:,7)/Fs>30)*100)
%     fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq)*100)
%     fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
%     fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
%     %     fprintf('... ... %g %% waves discarded because of pption neg elect\n',mean(all_Waves(:,16)>paramSW.min_pptionNeg)*100)
%     all_Waves(all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
%     %     all_Waves(all_Waves(:,16)>paramSW.min_pptionNeg | all_freq<paramSW.LimFrqW(1) | all_freq>paramSW.LimFrqW(2) | all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
%     
%     thr_Wave=[];
%     slow_Waves=[];
%     for nE=1:length(ChanLabels)
%         thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
%         temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
%         
%         thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
%         slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
%     end
    load([save_path filesep 'prct_DSS_SW_' SubID]);%,'slow_Waves','paramSW','ChanLabels')

    % ROI1 = importdata([ROIpath 'Preprobe10sARV_not_ot03422_s' SubID(end-2:end) '.mat']); % Bilat Ant cing
    % ROI2 = importdata([ROIpath 'Preprobe10sARV_not_ot2-3830_s' SubID(end-2:end) '.mat']); % Bilat mid-cing, L post cing
    % ROI3 = importdata([ROIpath 'Preprobe10sARV_not_ot_14-6430_s' SubID(end-2:end) '.mat']); %R cuenous, precuneous
    % ROI4 = importdata([ROIpath 'Preprobe10sARV_ot_not_36-42-4_s' SubID(end-2:end) '.mat']); % FFA

    % ROIfrontR = importdata([ROIpath 'Preprobe10sARV_front_mwmbVSot_roi_s' SubID(end-2:end) '.mat']); % marsbar
    % ROIfrontMP = importdata([ROIpath 'Preprobe10sARV_front_mwmbVSot_mpcorrected_s' SubID(end-2:end) '.mat']); % spm
    % ROIbackR = importdata([ROIpath 'Preprobe10sARV_back_mwmbVSot_roi_s' SubID(end-2:end) '.mat']); % marsmar
    % ROIbackMP = importdata([ROIpath 'Preprobe10sARV_back_mwmbVSot_mpcorrected_s' SubID(end-2:end) '.mat']); % spm

    % if exist([ROIpath 'Preprobe20sARV_front_mwmbVSot_spmcon4_s' SubID(end-2:end) '.mat'])==0
    %     continue;
    % end

     ROIfront = importdata([ROIpath 'Preprobe20sARV_front_mwmbVSot_mars_s' SubID(end-2:end) '.mat']); % marsbar  Preprobe20sARV_back_mwmbVSot_spmcon4
     % ROIfrontMP = importdata([ROIpath 'Preprobe20sARV_front_mwmbVSot_spmcon4filt_s' SubID(end-2:end) '.mat']); % filt
     ROIback = importdata([ROIpath 'Preprobe20sARV_back_mwmbVSot_mars_s' SubID(end-2:end) '.mat']); % marsmar
     % ROIbackMP = importdata([ROIpath 'Preprobe20sARV_back_mwmbVSot_spmcon4filt_s' SubID(end-2:end) '.mat']); % Preprobe5sARV_back_mwmbVSot_mars_s242


    


    for nP=unique(slow_Waves(:,2))'
        slow_Waves_perE=[];
        duration_of_probe=window_before_probes/60;
        
        temp_test_res=test_res(test_res(:,1)==probe_res(nP,5) & test_res(:,4)<=probe_res(nP,7),:);
        temp_go=temp_test_res(temp_test_res(:,5)~=3,:);
        temp_nogo=temp_test_res(temp_test_res(:,5)==3,:);

        if size(temp_go,1)<18 || size(temp_nogo,1)<2
            sprintf('... skipping probe %g in block %g from %s because missing data',probe_res(nP,1),probe_res(nP,5),SubID)
            continue;
        end
        gosize = 17; 
        nogosize = 1;

        temp_go=temp_go(end-gosize:end,:);
        temp_nogo=temp_nogo(end-nogosize:end,:);

        if temp_go(end,end)==0 % remove last miss 
            %sprintf('At probe %g in block %g from %s has a miss before probe',probe_res(nP,1),probe_res(nP,5),SubID)
            temp_go(end,:) =[];
        end
        
        for nE=1:length(ChanLabels)
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP)/duration_of_probe nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,5))/Fs)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,13)) NaN]];
        end
        table_length=size(table_SW,1);
        
        table_SW.SubID(table_length+(1:length(ChanLabels)))=repmat({SubID},length(ChanLabels),1);
        table_SW.Block(table_length+(1:length(ChanLabels)))=repmat(probe_res(nP,5),length(ChanLabels),1);
        table_SW.Probe(table_length+(1:length(ChanLabels)))=repmat(nP,length(ChanLabels),1);
        table_SW.Elec(table_length+(1:length(ChanLabels)))=ChanLabels;
        
        table_SW.SW_density(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,1);
        table_SW.SW_amplitude(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,2);
        table_SW.SW_frequency(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,3);
        table_SW.SW_downslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,4);
        table_SW.SW_upslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,5);
        

        Miss=mean(temp_go(:,end)==0);
        HitRT=nanmean(temp_go(:,11)-temp_go(:,8));
        FA=mean(temp_nogo(:,end-1)==0);
        
        this_State=States{probe_res(nP,17)};
        this_Vig=probe_res(nP,19);
        this_MWcat=probe_res(nP,18);
        
        table_SW.SART_Miss(table_length+(1:length(ChanLabels)))=repmat(Miss,length(ChanLabels),1);
        table_SW.SART_HitRT(table_length+(1:length(ChanLabels)))=repmat(HitRT,length(ChanLabels),1);
        table_SW.SART_FA(table_length+(1:length(ChanLabels)))=repmat(FA,length(ChanLabels),1);
        table_SW.SART_State(table_length+(1:length(ChanLabels)))=repmat({this_State},length(ChanLabels),1);
        table_SW.SART_Vig(table_length+(1:length(ChanLabels)))=repmat(this_Vig,length(ChanLabels),1);
        table_SW.SART_MWcat(table_length+(1:length(ChanLabels)))=repmat(this_MWcat,length(ChanLabels),1);

        probeOnset = probe_res(nP,2); 
        temp_test_res2=test_res(test_res(:,1)==probe_res(nP,5) & test_res(:,8)<=probeOnset,:);

        temp_go2=temp_test_res2(temp_test_res2(:,5)~=3,:);
        temp_nogo2=temp_test_res2(temp_test_res2(:,5)==3,:);

        temp_go2(:,14)=temp_go2(:,8)-probeOnset;
        temp_nogo2(:,14)=temp_nogo2(:,8)-probeOnset;

        temp_go2=temp_go2(temp_go2(:,14)>=-10,:);
        temp_nogo2=temp_nogo2(temp_nogo2(:,14)>=-10,:);

        if temp_go2(end,end)==0 % remove last miss 
            %sprintf('At probe %g in block %g from %s has a miss before probe',probe_res(nP,1),probe_res(nP,5),SubID)
            temp_go2(end,:) =[];
        end

        if strcmpi(this_State,{'ON'}) && (any(temp_nogo2(:,12)==0) || any(temp_go2(:,13)==0))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'ONerror'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'ON'}) && (all(temp_nogo2(:,12)==1) && all(temp_go2(:,13)==1))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'ONclean'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'MW'}) && (any(temp_nogo2(:,12)==0) || any(temp_go2(:,13)==0))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'MWerror'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'MW'}) && (all(temp_nogo2(:,12)==1) && all(temp_go2(:,13)==1))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'MWclean'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'MB'}) && (any(temp_nogo2(:,12)==0) || any(temp_go2(:,13)==0))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'MBerror'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'MB'}) && (all(temp_nogo2(:,12)==1) && all(temp_go2(:,13)==1))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'MBclean'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'DK'}) && (any(temp_nogo2(:,12)==0) || any(temp_go2(:,13)==0))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'DKerror'},length(ChanLabels),1);
        elseif strcmpi(this_State,{'DK'}) && (all(temp_nogo2(:,12)==1) && all(temp_go2(:,13)==1))
            table_SW.SART_StateError(table_length+(1:length(ChanLabels)))=repmat({'DKclean'},length(ChanLabels),1);
        end

        dataProbeind = unique(slow_Waves(slow_Waves(:,2)==nP,16));

        % this_ROI1=ROI1(dataProbeind,:);
        % table_SW.ROI1(table_length+(1:length(ChanLabels)))=repmat(this_ROI1,length(ChanLabels),1);
        % this_ROI2=ROI2(dataProbeind,:);
        % table_SW.ROI2(table_length+(1:length(ChanLabels)))=repmat(this_ROI2,length(ChanLabels),1);
        % this_ROI3=ROI3(dataProbeind,:);
        % table_SW.ROI3(table_length+(1:length(ChanLabels)))=repmat(this_ROI3,length(ChanLabels),1);
        % this_ROI4=ROI4(dataProbeind,:);
        % table_SW.ROI4(table_length+(1:length(ChanLabels)))=repmat(this_ROI4,length(ChanLabels),1);

        this_ROIfront=ROIfront(dataProbeind,:);
        table_SW.ROIfront(table_length+(1:length(ChanLabels)))=repmat(this_ROIfront,length(ChanLabels),1);
        % this_ROIfrontMP=ROIfrontMP(dataProbeind,:);
        % table_SW.ROIfrontMP(table_length+(1:length(ChanLabels)))=repmat(this_ROIfrontMP,length(ChanLabels),1);
        this_ROIback=ROIback(dataProbeind,:);
        table_SW.ROIback(table_length+(1:length(ChanLabels)))=repmat(this_ROIback,length(ChanLabels),1);
        % this_ROIbackMP=ROIbackMP(dataProbeind,:);
        % table_SW.ROIbackMP(table_length+(1:length(ChanLabels)))=repmat(this_ROIbackMP,length(ChanLabels),1);

        
    end
end
writetable(table_SW,[save_path filesep 'MW_MRI_SW_Behav_PerProbe_ROI_marsbar.txt']);


%%
addpath((path_fieldtrip))
ft_defaults;

% ChanLabels={EEG.chanlocs.labels};
ChanLabels(find(ismember(ChanLabels,'FPz')))={'Fpz'};
ChanLabels(find(ismember(ChanLabels,'FP1')))={'Fp1'};
table_SW.Elec(find(ismember(table_SW.Elec,'FPz')))={'Fpz'};
table_SW.Elec(find(ismember(table_SW.Elec,'FP1')))={'Fp1'};
table_SW.Elec=removecats(table_SW.Elec);

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

%% BEHAV ANALYSIS & PLOT for thesis 
% Correlation with Mindstate - ON vs OFF task 
table_SW.SART_ON=nan(size(table_SW,1),1);
table_SW.SART_ON(table_SW.SART_State=='ON')=1;
table_SW.SART_ON(table_SW.SART_State~='ON')=0;
table_SW.Block = double(table_SW.Block);
table_SW.Block = categorical(table_SW.Block, [1, 2, 3, 4], {'1st', '2nd', '3rd', '4th'}, 'Ordinal',true);

table_SW.SART_MW=nan(size(table_SW,1),1);
table_SW.SART_MW(table_SW.SART_State=='MW')=1;
table_SW.SART_MW(table_SW.SART_State=='ON')=0;

table_SW.SART_MB=nan(size(table_SW,1),1);
table_SW.SART_MB(table_SW.SART_State=='MB' | table_SW.SART_State=='DK')=1;
table_SW.SART_MB(table_SW.SART_State=='ON')=0;

table_SW.SART_Vig2=nan(size(table_SW,1),1);
table_SW.SART_Vig2(table_SW.SART_Vig==1 | table_SW.SART_Vig==2)=0;
table_SW.SART_Vig2(table_SW.SART_Vig==3 | table_SW.SART_Vig==4)=1;

table_SW.SART_ONerror=nan(size(table_SW,1),1);
table_SW.SART_ONerror(table_SW.SART_StateError=='ONerror')=1;
table_SW.SART_ONerror(table_SW.SART_StateError=='ONclean')=0;

table_SW.SART_MWerror=nan(size(table_SW,1),1);
table_SW.SART_MWerror(table_SW.SART_StateError=='MWerror')=1;
table_SW.SART_MWerror(table_SW.SART_StateError=='MWclean')=0;

table_SW.SART_MWinter=nan(size(table_SW,1),1);
table_SW.SART_MWinter(table_SW.SART_MWcat>1)=1; % personal or task
table_SW.SART_MWinter(table_SW.SART_MWcat==1)=0; % something in the room


subsettable_SW = table_SW(table_SW.Elec=={'Cz'},:); % to select 1 set of behavioral data

tabulate(subsettable_SW.SART_ON)
tabulate(subsettable_SW.SART_State)
tabulate(subsettable_SW.SART_MWinter)
tabulate(subsettable_SW.SART_StateError)

% Specify the variables of interest
groupVar = 'SART_State'; % Grouping variable
varsToAvg = {'SART_Miss', 'SART_FA', 'SART_HitRT'}; % Variables to calculate averages

% Use groupsummary to calculate averages per SART_ON
avgValues = groupsummary(subsettable_SW, groupVar, 'mean', varsToAvg);
avgValues.mean_SART_Miss=round(avgValues.mean_SART_Miss,2)*100;
avgValues.mean_SART_FA=round(avgValues.mean_SART_FA,2)*100;
avgValues.mean_SART_HitRT=round(avgValues.mean_SART_HitRT,3);
% Display the result
disp(avgValues);

% Use groupsummary to calculate averages per SART_ON
avgValues = groupsummary(subsettable_SW, groupVar, 'std', varsToAvg);
avgValues.std_SART_Miss=round(avgValues.std_SART_Miss,2)*100;
avgValues.std_SART_FA=round(avgValues.std_SART_FA,2)*100;
avgValues.std_SART_HitRT=round(avgValues.std_SART_HitRT,3);
% Display the result
disp(avgValues);


lme_miss = fitlme(subsettable_SW, 'SART_Miss ~ 1 + Block + SART_ON + (1+Block|SubID)');
lme_rt = fitlme(subsettable_SW, 'SART_HitRT ~ 1 + Block + SART_ON + (1+Block|SubID)');
lme_fa = fitlme(subsettable_SW, 'SART_FA ~ 1 + Block + SART_ON + (1+Block|SubID)');


lme_miss = fitglme(subsettable_SW, 'SART_ON ~ 1 + Block + SART_Miss + (1+Block|SubID)' ,'Distribution','binomial');
lme_rt = fitglme(subsettable_SW, ' SART_ON ~ 1 + Block + SART_HitRT + (1+Block|SubID)' ,'Distribution','binomial');
lme_fa = fitglme(subsettable_SW, ' SART_ON ~ 1 + Block + SART_FA + (1+Block|SubID)','Distribution','binomial');


lme_miss = fitlme(subsettable_SW, 'SART_Miss ~ 1 + Block + SART_MWinter + (1+Block|SubID)');
lme_rt = fitlme(subsettable_SW, 'SART_HitRT ~ 1 + Block + SART_MWinter + (1+Block|SubID)');
lme_fa = fitlme(subsettable_SW, 'SART_FA ~ 1 + Block + SART_MWinter + (1+Block|SubID)');


lme_miss = fitlme(subsettable_SW, 'SART_Miss ~ 1 + Block + SART_MWerror + (1+Block|SubID)');
lme_rt = fitlme(subsettable_SW, 'SART_HitRT ~ 1 + Block + SART_MWerror + (1+Block|SubID)');
lme_fa = fitlme(subsettable_SW, 'SART_FA ~ 1 + Block + SART_MWerror + (1+Block|SubID)');



subsettable = subsettable_SW(:,{'SubID','Block','Probe','SART_Miss','SART_FA','SART_HitRT','SART_ON','SART_State','SART_Vig','SART_MWcat','SART_StateError'});

writetable(subsettable,[save_path filesep 'MW_MRI_ProbeOnly.txt']);


% lme_numeric_sl = fitlme(table_SW, 'SW_density ~ NumBlock + (1 + NumBlock|SubID)');
% % Ordinal model
% lme_ordinal = fitlme(table_SW, 'SW_density ~ Block + (1|SubID)');
% lme_ordinal_sl = fitlme(table_SW, 'SW_density ~ Block + (1+Block|SubID)');
% 
% comparison_result = compare(lme_numeric, lme_ordinal, lme_numeric_sl, lme_ordinal_sl);
% comparison_result = compare(lme_numeric, lme_numeric_sl);
% comparison_result = compare(lme_ordinal, lme_ordinal_sl);
% disp(comparison_result);


%% Overall correlation between SW distibution and performance/mindstate  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FIG 5
% Correlation with behaviour
clear *_effect
newlabels=layout.label(1:end-2);

% Correlation with Mindstate - ON vs OFF task 
table_SW.SART_ON=nan(size(table_SW,1),1);
table_SW.SART_ON(table_SW.SART_State=='ON')=1;
table_SW.SART_ON(table_SW.SART_State~='ON')=0;
table_SW.Block = double(table_SW.Block);
table_SW.Block = categorical(table_SW.Block, [1, 2, 3, 4], {'1st', '2nd', '3rd', '4th'}, 'Ordinal',true);



for nCh=1:length(newlabels)
    mdl_ON_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+Block+SW_density+(1+Block|SubID)','Distribution','binomial');
    ON_D_effect(nCh,1)=mdl_ON_D.Coefficients.tStat(match_str(mdl_ON_D.CoefficientNames,'SW_density'));
    ON_D_effect(nCh,2)=mdl_ON_D.Coefficients.pValue(match_str(mdl_ON_D.CoefficientNames,'SW_density'));

    mdl_ON_A=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+Block+SW_amplitude+(1+Block|SubID)','Distribution','binomial');
    ON_A_effect(nCh,1)=mdl_ON_A.Coefficients.tStat(match_str(mdl_ON_A.CoefficientNames,'SW_amplitude'));
    ON_A_effect(nCh,2)=mdl_ON_A.Coefficients.pValue(match_str(mdl_ON_A.CoefficientNames,'SW_amplitude'));

    mdl_ON_F=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+Block+SW_frequency+(1+Block|SubID)','Distribution','binomial');
    ON_F_effect(nCh,1)=mdl_ON_F.Coefficients.tStat(match_str(mdl_ON_F.CoefficientNames,'SW_frequency'));
    ON_F_effect(nCh,2)=mdl_ON_F.Coefficients.pValue(match_str(mdl_ON_F.CoefficientNames,'SW_frequency'));

    mdl_ON_Do=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+Block+SW_downslope+(1+Block|SubID)','Distribution','binomial');
    ON_Do_effect(nCh,1)=mdl_ON_Do.Coefficients.tStat(match_str(mdl_ON_Do.CoefficientNames,'SW_downslope'));
    ON_Do_effect(nCh,2)=mdl_ON_Do.Coefficients.pValue(match_str(mdl_ON_Do.CoefficientNames,'SW_downslope'));

    mdl_ON_U=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+Block+SW_upslope+(1+Block|SubID)','Distribution','binomial');
    ON_U_effect(nCh,1)=mdl_ON_U.Coefficients.tStat(match_str(mdl_ON_U.CoefficientNames,'SW_upslope'));
    ON_U_effect(nCh,2)=mdl_ON_U.Coefficients.pValue(match_str(mdl_ON_U.CoefficientNames,'SW_upslope'));
    fprintf('Electrode %s\n',newlabels{nCh})

end



% for the thesis
table_SW.Block = double(table_SW.Block);

for nCh=1:length(newlabels)
    mdl_ON_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+SW_density+(1+Block|SubID)','Distribution','binomial');
    ON_D_effect(nCh,1)=mdl_ON_D.Coefficients.tStat(match_str(mdl_ON_D.CoefficientNames,'SW_density'));
    ON_D_effect(nCh,2)=mdl_ON_D.Coefficients.pValue(match_str(mdl_ON_D.CoefficientNames,'SW_density'));
    fprintf('Electrode %s\n',newlabels{nCh})
    
    % mdl_ON_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+SW_amplitude+(1|SubID)','Distribution','binomial'); % df 970
    % ON_D_effect(nCh,1)=mdl_ON_D.Coefficients.tStat(match_str(mdl_ON_D.CoefficientNames,'SW_amplitude'));
    % ON_D_effect(nCh,2)=mdl_ON_D.Coefficients.pValue(match_str(mdl_ON_D.CoefficientNames,'SW_amplitude'));
    % fprintf('Electrode %s\n',newlabels{nCh})

end



% % Numeric model
% lme_numeric = fitlme(table_SW, 'SW_density ~ NumBlock + (1|SubID)');
% lme_numeric_sl = fitlme(table_SW, 'SW_density ~ NumBlock + (1 + NumBlock|SubID)');
% % Ordinal model
% lme_ordinal = fitlme(table_SW, 'SW_density ~ Block + (1|SubID)');
% lme_ordinal_sl = fitlme(table_SW, 'SW_density ~ Block + (1+Block|SubID)');
% 
% comparison_result = compare(lme_numeric, lme_ordinal, lme_numeric_sl, lme_ordinal_sl);
% comparison_result = compare(lme_numeric, lme_numeric_sl);
% comparison_result = compare(lme_ordinal, lme_ordinal_sl);
% disp(comparison_result);


% lme_full = fitglme(table_SW, 'SART_ON ~ 1 + Block + SW_density + (1 + Block | SubID) + (1 | Elec)', ...
%                    'Distribution', 'binomial');
% lme_reduced = fitglme(table_SW, 'SART_ON ~ 1 + SW_density + (1 + Block | SubID) + (1 | Elec)', ...
%                       'Distribution', 'binomial');
% comparison_result = compare(lme_reduced, lme_full);
% disp(comparison_result);



%% Figure

effect = [ON_D_effect(:,1), ON_A_effect(:,1), ON_F_effect(:,1), ON_Do_effect(:,1), ON_U_effect(:,1)];
all_pV = [ON_D_effect(:,2); ON_A_effect(:,2); ON_F_effect(:,2); ON_Do_effect(:,2); ON_U_effect(:,2)];
%all_pV = [ON_D_effect(:,2)];

% effect = [ROI1_effect(:,1), ROI1_effect(:,3), ROI1_effect(:,5)];
% all_pV = [ROI1_effect(:,2), ROI1_effect(:,4), ROI1_effect(:,6)];
FDR_Thr=0.05;%fdr(all_pV,0.05);
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
minmax = ceil(max(max(abs(effect))));


figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Set the figure position and size
simpleTopoPlot_ft(ON_D_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_D_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
%ft_plot_lay_me(layout, 'chanindx', find(abs(ON_D_effect(:,1))>1.7109), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON Density block intercept', 'FontSize', 16)

%
figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Set the figure position and size
subplot(2,2,1)
simpleTopoPlot_ft(ON_D_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_D_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON Density', 'FontSize', 16)

subplot(2,2,2);
simpleTopoPlot_ft(ON_A_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_A_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON Amplitude', 'FontSize', 16)

subplot(2,2,3);
simpleTopoPlot_ft(ON_Do_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_Do_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON Downslope', 'FontSize', 16)

subplot(2,2,4);
simpleTopoPlot_ft(ON_U_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_U_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON Upslope', 'FontSize', 16)

%% ROI
newlabels=layout.label(1:end-2);

%table_SW.Block = categorical(table_SW.Block, [1, 2, 3, 4], {'1st', '2nd', '3rd', '4th'}, 'Ordinal',true);
%table_SW.Block = double(table_SW.Block);
table_SW.Block = categorical(table_SW.Block);



% ROIfrontR ROIfrontMP ROIbackR ROIbackMP ROI3 ROI14

clear ROI_effect ROI_A_effect
for nCh=1:length(newlabels)
    mdl_ROI1_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROIbackR~1+SW_density*Block + (1+Block|SubID)','Distribution','Normal');
    ROI_effect(nCh,1)=mdl_ROI1_D.Coefficients.tStat(match_str(mdl_ROI1_D.CoefficientNames,'SW_density')); % SW_density Block Block:SW_density
    ROI_effect(nCh,2)=mdl_ROI1_D.Coefficients.pValue(match_str(mdl_ROI1_D.CoefficientNames,'SW_density')); % Block_2nd Block_2nd:SW_density (Intercept) Block_3rd Block_4th
    fprintf('Electrode %s\n',newlabels{nCh})

    mdl_ROI1_A=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROIbackR~1+SW_amplitude*Block + (1+Block|SubID)','Distribution','Normal');
    ROI_A_effect(nCh,1)=mdl_ROI1_A.Coefficients.tStat(match_str(mdl_ROI1_A.CoefficientNames,'SW_amplitude'));
    ROI_A_effect(nCh,2)=mdl_ROI1_A.Coefficients.pValue(match_str(mdl_ROI1_A.CoefficientNames,'SW_amplitude'));

end
  

effect = [ROI_effect(:,1), ROI_A_effect(:,1)];
FDR_Thr=0.05; %fdr(all_pV,0.05);
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
minmax = ceil(max(max(abs(effect))));

figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Set the figure position and size
subplot(1,2,1)
simpleTopoPlot_ft(ROI_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ROI_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI:Back - Density (SPM CON  -20 sec- Rand eff. Block=cat. Fixed with int) ', 'FontSize', 16) %  (Marsbar - bp filtered  - 5 sec - Block=continous)

subplot(1,2,2)
simpleTopoPlot_ft(ROI_A_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ROI_A_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI:Back - Amplitude (SPM CON -20 sec- Rand eff. Block=cat. Fixed with int) ', 'FontSize', 16)

%%

figure;
hold on;
for nCh=1:length(newlabels)
    mdl_ROI1_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROIbackR~1+SW_density+(1+Block|SubID)','Distribution','Normal');
    fprintf('Electrode %s\n',newlabels{nCh})


% Extract fixed effect coefficients
intercept = mdl_ROI1_D.Coefficients.Estimate(1); % Intercept
slope = mdl_ROI1_D.Coefficients.Estimate(2);     % Slope for SW_density

% Generate a range of SW_density values
sw_density_values = linspace(min(table_SW.SW_density), max(table_SW.SW_density), 100);

% Predict ROIbackR using the fixed effects
predicted_ROIbackR = intercept + slope * sw_density_values;

% Plot the effect of SW_density
plot(sw_density_values, predicted_ROIbackR, 'b-', 'LineWidth', 2);
xlabel('SW_density');
ylabel('Predicted ROIbackR');
title('Effect of SW_density on ROIbackR');
ylim([150 154])
grid on;

% Add confidence intervals (Optional)

ci_upper = predicted_ROIbackR + 1.96 * mdl_ROI1_D.Coefficients.SE(2); % Upper CI
ci_lower = predicted_ROIbackR - 1.96 * mdl_ROI1_D.Coefficients.SE(2); % Lower CI
fill([sw_density_values, fliplr(sw_density_values)], ...
     [ci_upper, fliplr(ci_lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
legend('Predicted ROIbackR', '95% Confidence Interval');
%hold off;
end



figure;
% Create a subplot
subplot(2, 2, 4); % 2x2 grid of subplots
hold on;
for nCh=1:length(newlabels)
    mdl_ROI1_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROIbackR~1+SW_density+(1+Block|SubID)','Distribution','Normal');
    fprintf('Electrode %s\n',newlabels{nCh})
% Extract fixed effect coefficients
intercept = mdl_ROI1_D.Coefficients.Estimate(1); % Intercept
slope = mdl_ROI1_D.Coefficients.Estimate(2);     % Slope for SW_density

% Unique block levels
block_levels = [1, 2, 3, 4]; % Numeric representation of blocks

% Extract random effects
random_effects = randomEffects(mdl_ROI1_D); % 1x100 vector
random_effects = reshape(random_effects, [], 2); % Reshape to match SubID (Intercept and Block slopes)

% Generate a range of SW_density values
sw_density_values = linspace(min(table_SW.SW_density), max(table_SW.SW_density), 100);

% Loop through each block and create a subplot
for i = 1:length(block_levels)
    % Random effect intercept for the current block
    random_intercept_adjustment = mean(random_effects(:, 1)); % Mean intercept adjustment across subjects
    random_block_slope = mean(random_effects(:, 2));          % Mean slope for the block
    
    % Compute adjusted intercept for the current block
    adjusted_intercept = intercept + random_intercept_adjustment + random_block_slope * block_levels(i);
    
    % Predict ROIbackR
    predicted_ROIbackR = adjusted_intercept + slope * sw_density_values;
    
    subplot(2, 2, i); % 2x2 grid of subplots
    hold on
    plot(sw_density_values, predicted_ROIbackR, 'b-', 'LineWidth', 2);
    xlabel('SW_density');
    ylabel('Predicted ROIbackR');
    title(sprintf('Block %d', block_levels(i)));
    ylim([150 152])
    grid on;
end

% Add a super title for the entire figure
sgtitle('Effect of SW_density on ROIbackR Across Blocks');

end












% figure;
% simpleTopoPlot_ft(ON_effect(:,1), layout,'on',[],0,1);
% ft_plot_lay_me(layout, 'chanindx', find(ON_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
% colorbar;
% colormap(cmap2);
% caxis([-1 1]*5)
% title('OFF vs ON', 'FontSize', 16)
% c = colorbar; c.Label.String = 't-value'; c.Label.FontSize = 14; c.Label.Rotation = 270; c.Label.Position(1) = 4; c.FontSize = 14;
% 


% %% Topography
% % ChanLabels={EEG.chanlocs.labels};
% ChanLabels(find(ismember(ChanLabels,'FPz')))={'Fpz'};
% ChanLabels(find(ismember(ChanLabels,'FP1')))={'Fp1'};
% 
% cfg = [];
% cfg.layout = 'elec1005.lay';
% cfg.center      = 'yes';
% cfg.channel = ChanLabels(~ismember(ChanLabels,{'TP9','TP10'}));
% layout=ft_prepare_layout(cfg);
% % layout.label(match_str(layout.label,{'TP9','TP10'}))=[];
% correspCh=[];
% for nCh=1:length(layout.label)-2
%     correspCh(nCh)=match_str(ChanLabels,layout.label{nCh});
% end
% 
% figure;
% subplot(2,2,1);
% topo_plot=squeeze(nanmean(nanmean(SW_dens_perProbe,2),1));
% % topo_plot=squeeze(nanmean(SW_dens,1));
% % topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=0;
% simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
% colorbar;
% title('SW density')
% 
% % figure;
% subplot(2,2,2);
% topo_plot=squeeze(nanmean(all_P2P,1));
% topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
% simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
% colorbar;
% title('SW amplitude')
% 
% subplot(2,2,3);
% topo_plot=squeeze(nanmean(all_PosSl,1));
% topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
% simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
% colorbar;
% title('SW positive slope')
% 
% subplot(2,2,4);
% topo_plot=squeeze(nanmean(all_NegSl,1));
% topo_plot(match_str(ChanLabels,{'TP9','TP10'}))=NaN;
% simpleTopoPlot_ft(topo_plot(correspCh), layout,'labels',[],0,1);
% colorbar;
% title('SW negative slope')

%%

%  Fit a GLME 
lme = fitglme(table_SW, 'SART_ON ~ Block + (1|SubID) + (1 | Elec)', 'Distribution', 'binomial');
disp('Mixed-Effects Model Summary:');
disp(lme);

% Extract and interpret Block effects from GLME
fixed_effects = lme.Coefficients; % Extract fixed effects
disp('Fixed Effects from GLME:');
disp(fixed_effects);

% Convert log-odds to probabilities
log_odds = fixed_effects.Estimate; % Fixed effects for Block
odds = exp(log_odds); % Convert to odds
probs = odds ./ (1 + odds); % Convert to probabilities
fprintf('Predicted probabilities:\n');
disp(probs);

% Visualize predicted probabilities by Block
predicted_probs = fitted(lme); % Predicted probabilities
figure;
boxchart(table_SW.Block, predicted_probs, 'MarkerStyle', 'o', 'Notch', 'on');
xlabel('Block');
ylabel('Predicted Probability of SART_ON');
title('Effect of Block on SART_ON');
grid on;


%%

% Histogram with Normal Distribution Overlay
figure;
histogram(table_SW.ROIfrontR, 'Normalization', 'pdf');
hold on;
x = linspace(min(table_SW.ROIfrontR), max(table_SW.ROIfrontR), 100);
y = normpdf(x, mean(table_SW.ROIfrontR, 'omitnan'), std(table_SW.ROIfrontR, 'omitnan'));
plot(x, y, 'r', 'LineWidth', 2);
title('Histogram with Normal Distribution Overlay');
xlabel('ROIfrontR');
ylabel('Density');
hold off;

% Histogram with Normal Distribution Overlay
figure;
histogram(table_SW.ROIbackR, 'Normalization', 'pdf');
hold on;
x = linspace(min(table_SW.ROIbackR), max(table_SW.ROIbackR), 100);
y = normpdf(x, mean(table_SW.ROIbackR, 'omitnan'), std(table_SW.ROIbackR, 'omitnan'));
plot(x, y, 'r', 'LineWidth', 2);
title('Histogram with Normal Distribution Overlay');
xlabel('ROIbackR');
ylabel('Density');

%%

parpool;
tic
% Identify unique combinations of SubID, Block, and Probe
[groupIndices, uniqueGroups] = findgroups(table_SW(:, {'SubID', 'Block', 'Probe'}));
% Pre-allocate output variables
numPermutations = 1000;
numChannels = length(newlabels);
ON_D_Perm_tStat = nan(numChannels, numPermutations);
ON_D_Perm_pVal = nan(numChannels, numPermutations);

parfor nperm=1:numPermutations
tempTable = table_SW; % Create a temporary copy for parallel execution

tempTable.Block = double(tempTable.Block);
% ROI
tempTable.ROIfrontRPERM = nan(height(table_SW), 1);
%get all ROI values
orginalVals = unique(tempTable.ROIfrontR);
%orginalVals = unique(tempTable.ROIbackR);

orginalValsnoreuse = orginalVals;
% Loop through each group (unique Probe + Elec) to shuffle ROIfrontR
    for g = 1:size(orginalVals, 1)
        probeRows = groupIndices == g;
        permutedValue = orginalValsnoreuse(randperm(length(orginalValsnoreuse), 1)); % Pick one random value
        % Assign the permuted value to all rows in this group
        tempTable.ROIfrontRPERM(probeRows) = permutedValue;
        orginalValsnoreuse(orginalValsnoreuse(:,1)==permutedValue)=[];
    end

% % Correlation with Mindstate - ON vs OFF task 
% tempTable.SART_ONPERM=nan(size(tempTable,1),1);
% % Generate random 0 or 1 for each unique group
% randomValues = randi([0, 1], size(uniqueGroups, 1), 1);
% % Map the random values back to the table
% tempTable.SART_ONPERM = randomValues(groupIndices);

% this section takes a lot of time -- Add later
% block_effect = [0.6403,0.3338,0.3336,0.4702]; %  log-odds for Block levels 1st, 2nd, 3rd, 4th
% % Convert categorical directly to numeric
% numeric_blocks = double(tempTable.Block);
% % Generate probabilities based on Block effect
% logit_probs = block_effect(numeric_blocks); % Get log-odds for each row based on Block
% probssim = 1 ./ (1 + exp(-logit_probs)); % Convert log-odds to probabilities
% % Simulate SART_ONPERM based on probabilities
% tempTable.SART_ONPERM = rand(size(tempTable, 1), 1) <= probssim';

% Temporary variables to store results for each permutation
tStatTemp = nan(numChannels, 1);
pValTemp = nan(numChannels, 1);

    for nCh=1:length(newlabels)   
        mdl_ON_D_perm=fitglme(tempTable(tempTable.Elec==newlabels{nCh},:),'ROIfrontRPERM~1+SW_density+(1+Block|SubID)','Distribution','Normal');
        %mdl_ON_D_perm=fitglme(tempTable(tempTable.Elec==newlabels{nCh},:),'SART_ONPERM~1+SW_amplitude+(1|SubID)','Distribution','binomial');
        %mdl_ON_D_perm=fitglme(tempTable(tempTable.Elec==newlabels{nCh},:),'SART_ONPERM~1+SW_density+(1|SubID)','Distribution','binomial');
        %mdl_ON_D_perm=fitglme(tempTable(tempTable.Elec==newlabels{nCh},:),'SART_ONPERM~1+Block+SW_density+(1+Block|SubID)','Distribution','binomial');
        % tStatTemp(nCh)=mdl_ON_D_perm.Coefficients.tStat(match_str(mdl_ON_D_perm.CoefficientNames,'SW_amplitude'));
        % pValTemp(nCh)=mdl_ON_D_perm.Coefficients.pValue(match_str(mdl_ON_D_perm.CoefficientNames,'SW_amplitude'));
        tStatTemp(nCh)=mdl_ON_D_perm.Coefficients.tStat(match_str(mdl_ON_D_perm.CoefficientNames,'SW_density'));
        pValTemp(nCh)=mdl_ON_D_perm.Coefficients.pValue(match_str(mdl_ON_D_perm.CoefficientNames,'SW_density'));
        %fprintf('Electrode %s\n',newlabels{nCh})
    end
    % Store temporary results into the main matrices
    ON_D_Perm_tStat(:, nperm) = tStatTemp;
    ON_D_Perm_pVal(:, nperm) = pValTemp;
    fprintf('Running permutation %d\n',nperm)
end
toc
% Close the parallel pool
delete(gcp);

%%
% CHECK SIM

% without parfor run time 1.1 min per permutation
% with parfor 0.3 min per permutation
% with parfor and simple model without block: 8.5812 min for 1000
% permutation
% i remove random block slope

figure;
histogram(reshape(ON_D_Perm_tStat, 1, []))

% Flatten the p-values matrix into a single vector for plotting
flat_p_values = reshape(ON_D_Perm_pVal, 1, []);

% Plot the histogram of p-values
figure;
histogram(flat_p_values, 'BinWidth', 0.01);
hold on;

% Add a vertical line at p = 0.05
xline(0.05, 'r--', 'LineWidth', 2, 'DisplayName', 'p = 0.05');

% Calculate the percentage of p-values below 0.05
percent_below_005 = mean(flat_p_values < 0.05) * 100;

% Add the percentage as a subtitle
subtitle(sprintf('%.2f%% of p-values are below 0.05', percent_below_005), 'FontSize', 12, 'Color', 'red');


% Customize plot
xlabel('p-value');
ylabel('Probability');
title('Histogram of Permutation p-values');
legend('Histogram', 'p = 0.05');
grid on;
hold off;


% check an example 
thr=tinv(0.975, 25 - 1); % Two-tailed t-threshold
figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Set the figure position and size
simpleTopoPlot_ft(ON_D_Perm_tStat(:,2), layout,'on',[],0,1);
%ft_plot_lay_me(layout, 'chanindx', find(ON_D_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
ft_plot_lay_me(layout, 'chanindx', find(abs(ON_D_Perm_tStat(:,2))>thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON Density', 'FontSize', 16)


%save('ROIfront_nperm100_density_withblock_null.mat', 'ON_D_Perm_tStat', 'ON_D_Perm_pVal');



%%

% Permutation
% we will determine the channel neighbours
cfg = [];
cfg.layout = 'elec1005.lay';
cfg.center      = 'yes';
cfg.channel = ChanLabels(~ismember(ChanLabels,{'TP9','TP10'}));
layout=ft_prepare_layout(cfg);
cfg.method        = 'triangulation';
nb = ft_prepare_neighbours(cfg);

electrode_neighbors=nb.neighblabel;
electrode_label=nb.label;


% Create a neighbors lookup table
electrodes = {nb.label}; % Extract electrode labels
neighbors = cell(size(nb));

for i = 1:length(nb)
    % Find the indices of the neighboring electrodes
    neighbors{i} = find(ismember(electrodes, nb(i).neighblabel));
end

t_values = ON_D_effect(:,1);

% Define threshold for cluster formation (e.g., corresponding to p = 0.05)
threshold = tinv(0.975, 25 - 1); % Two-tailed t-threshold
% Identify suprathreshold electrodes
suprathreshold = abs(t_values) > threshold;

% Initialize cluster labels
cluster_labels = zeros(size(t_values));
current_label = 0;

for e = 1:length(t_values)
    if suprathreshold(e) && cluster_labels(e) == 0
        % Start a new cluster
        current_label = current_label + 1;
        stack = [e]; % Stack for flood-fill algorithm
        
        while ~isempty(stack)
            current = stack(end);
            stack(end) = []; % Pop the stack
            
            if cluster_labels(current) == 0 && suprathreshold(current)
                % Label the current electrode
                cluster_labels(current) = current_label;
                
                % Add neighbors to the stack
                for neighbor = neighbors{current}
                    if cluster_labels(neighbor) == 0 && suprathreshold(neighbor)
                        stack = [stack, neighbor];
                    end
                end
            end
        end
    end
end



unique_clusters = unique(cluster_labels);
unique_clusters(unique_clusters == 0) = []; % Exclude background

% Compute cluster masses
cluster_masses = arrayfun(@(c) sum(t_values(cluster_labels == c)), unique_clusters);


for p = 1:numPermutations
% Identify clusters in permuted data
    perm_t_values = ON_D_Perm_tStat(:,p);
    perm_suprathreshold = abs(perm_t_values) > threshold;
    perm_cluster_labels = zeros(size(perm_t_values));
    current_label = 0;

    for e = 1:length(perm_t_values)
        if perm_suprathreshold(e) && perm_cluster_labels(e) == 0
            current_label = current_label + 1;
            stack = [e];
            
            while ~isempty(stack)
                current = stack(end);
                stack(end) = [];
                
                if perm_cluster_labels(current) == 0 && perm_suprathreshold(current)
                    perm_cluster_labels(current) = current_label;
                    
                    for neighbor = neighbors{current}
                        if perm_cluster_labels(neighbor) == 0 && perm_suprathreshold(neighbor)
                            stack = [stack, neighbor];
                        end
                    end
                end
            end
        end
    end

% Compute cluster mass for permuted data
    unique_perm_clusters = unique(perm_cluster_labels);
    unique_perm_clusters(unique_perm_clusters == 0) = [];
    perm_cluster_masses = arrayfun(@(c) sum(perm_t_values(perm_cluster_labels == c)), unique_perm_clusters);

    % Store the largest cluster mass
    if ~isempty(perm_cluster_masses)
        null_cluster_masses(p) = max(perm_cluster_masses);
    end

end

% Compute p-values for observed clusters
p_values = arrayfun(@(mass) mean(abs(null_cluster_masses) >= mass), abs(cluster_masses));

null_cluster_masses(abs(null_cluster_masses)>=abs(cluster_masses(2,:)),:)

% Display significant clusters
significant_clusters = find(p_values < 0.05);
fprintf('Significant clusters:\n');
for c = significant_clusters
    fprintf('Cluster %d: p = %.3f\n', c, p_values(c));
end


%%
% https://benediktehinger.de/blog/science/threshold-free-cluster-enhancement-explained/

% we will determine the channel neighbours
cfg = [];
cfg.layout = 'elec1005.lay';
cfg.center      = 'yes';
cfg.channel = ChanLabels(~ismember(ChanLabels,{'TP9','TP10'}));
layout=ft_prepare_layout(cfg);
cfg.method        = 'triangulation';
nb = ft_prepare_neighbours(cfg);

electrode_neighbors=nb.neighblabel;
electrode_label=nb.label;


%%
% Create a neighbors lookup table
electrodes = {nb.label}; % Extract electrode labels
neighbors = cell(size(nb));

for i = 1:length(nb)
    % Find the indices of the neighboring electrodes
    neighbors{i} = find(ismember(electrodes, nb(i).neighblabel));
end


% Parameters for TFCE
%t_values = ON_D_effect(:,1);
t_values = ROI_effect(:,1);
%t_values = ROI_A_effect(:,1);

H = 2; % Height weight -  2
E = 0.5; % Extent weight -  0.5
nThresholds = 10; % Number of thresholds to evaluate
tMax = max(abs(t_values)); % Maximum t-value for thresholds
thresholds = linspace(0, tMax, nThresholds); % Define thresholds

% Initialize TFCE values
%tfce_values = zeros(size(thresholds,2), length(t_values));
tfce_values = zeros(1, length(t_values));
ind=1;

% Compute TFCE for each electrode
for h = thresholds
    % Identify suprathreshold electrodes
    suprathreshold = abs(t_values) > h;

    % Label connected clusters
    cluster_labels = zeros(size(t_values));
    current_label = 0;
    for e = 1:length(t_values)
        if suprathreshold(e) && cluster_labels(e) == 0
            % Start a new cluster
            current_label = current_label + 1;
            stack = [e]; % Flood-fill algorithm

            while ~isempty(stack)
                current = stack(end);
                stack(end) = []; % Pop the stack
                if cluster_labels(current) == 0 && suprathreshold(current)
                    cluster_labels(current) = current_label;

                    % Add neighbors to the stack
                    for neighbor = neighbors{current}
                        if cluster_labels(neighbor) == 0 && suprathreshold(neighbor)
                            stack = [stack, neighbor];
                        end
                    end
                end
            end
        end
    end

    % Compute extent for each cluster
    unique_clusters = unique(cluster_labels);
    unique_clusters(unique_clusters == 0) = []; % Exclude background


    for cluster = unique_clusters'

        if isempty(cluster)
            continue; % Skip if the electrode has no neighbors
        end
        cluster_indices = find(cluster_labels == cluster);
        cluster_extent = numel(cluster_indices);

        % Enhance TFCE values for electrodes in the cluster
        % tfce_values(ind,cluster_indices) = tfce_values(ind,cluster_indices) + ...
        %     (h^H) * (cluster_extent^E) * (thresholds(2) - thresholds(1));
        tfce_values(cluster_indices) = tfce_values(cluster_indices) + ...
            (h^H) * (cluster_extent^E) * (thresholds(2) - thresholds(1)); %evenly spaced!
    end
    %ind=ind+1;
end

null_tfce_values = zeros(length(t_values), numPermutations);
for p = 1:numPermutations
% Identify clusters in permuted data
    perm_t_values = ON_D_Perm_tStat(:,p);
    %ind=1;
    % Compute TFCE for permuted data
    perm_tfce = zeros(1, length(perm_t_values));

    for h = thresholds
        perm_suprathreshold = abs(perm_t_values) > h;
        perm_cluster_labels = zeros(size(perm_t_values));
        current_label = 0;

        for e = 1:length(perm_t_values)
            if perm_suprathreshold(e) && perm_cluster_labels(e) == 0
                current_label = current_label + 1;
                stack = [e];
                
                while ~isempty(stack)
                    current = stack(end);
                    stack(end) = [];
                    
                    if perm_cluster_labels(current) == 0 && perm_suprathreshold(current)
                        perm_cluster_labels(current) = current_label;
                        
                        for neighbor = neighbors{current}
                            if perm_cluster_labels(neighbor) == 0 && perm_suprathreshold(neighbor)
                                stack = [stack, neighbor];
                            end
                        end
                    end
                end
            end
        end
 

   % Compute extent and enhance permuted TFCE values
        unique_perm_clusters = unique(perm_cluster_labels);
        unique_perm_clusters(unique_perm_clusters == 0) = [];
        for cluster = unique_perm_clusters'
            if isempty(cluster)
                continue; % Skip if cluster empty
            end
            cluster_indices = find(perm_cluster_labels == cluster);
            cluster_extent = numel(cluster_indices);
            % perm_tfce(ind,cluster_indices) = perm_tfce(ind,cluster_indices) + ...
            %     (h^H) * (cluster_extent^E) * (thresholds(2) - thresholds(1));
            perm_tfce(cluster_indices) = perm_tfce(cluster_indices) + ...
              (h^H) * (cluster_extent^E) * (thresholds(2) - thresholds(1));
        end
        %=ind+1;
    end
  % Store TFCE values for this permutation
    null_tfce_values(:, p) = perm_tfce;
end


%%
% Check
% Compute p-values
p_values = mean(null_tfce_values >= tfce_values', 2);

% Compute p-values for observed TFCE values
%p_values = mean(null_tfce_values(:,1) >= tfce_values(7,:)', 2);

% Display significant electrodes
significant_electrodes = find(p_values < 0.05);

fprintf('Significant electrodes:\n');
for i = 1:length(significant_electrodes)
    e = significant_electrodes(i); % Electrode index
    fprintf('Electrode %d: p = %.3f \n', e, p_values(e));
end



figure;
set(gcf, 'Position', [100, 100, 500, 500]); % Set the figure position and size
simpleTopoPlot_ft( ROI_effect(:,1), layout,'on',[],0,1);
%simpleTopoPlot_ft( ON_D_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', significant_electrodes, 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI:Front Density', 'FontSize', 16)
format_fig;


