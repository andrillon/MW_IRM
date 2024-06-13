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
files=dir([data_path filesep filesep 'MWMRI*clean5.set']);


%%
table_SW=array2table(zeros(0,15),'VariableNames',{'SubID','Block','Probe','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SART_Miss','SART_FA','SART_HitRT','SART_ON','SART_State','SART_Vig'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.Elec=categorical(table_SW.Elec);
table_SW.SART_State=categorical(table_SW.SART_State);

States={'ON','MW','MB','DK'};
%% loop across trials for baseline blocks
nFc=0;
all_SW_probes=[];
window_before_probes=20; % in seconds
all_table_behav_SW=[];
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
    if nFc==1
        addpath(genpath(path_eeglab));
        EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
        rmpath(genpath(path_eeglab));
    end

    Fs=500;
  load([save_path filesep 'prct_DSS_SW_delta_' SubID])

    table_behav_SW=array2table(zeros(0,10+length(ChanLabels)),'VariableNames',[{'Block','Probe','Trial','Condition','DistToProbe','RT','FA','Miss','MS','Vigilance'},ChanLabels]);
    for nP=unique(slow_Waves(:,2))'
        temp_test_res=test_res(test_res(:,1)==probe_res(nP,5) & test_res(:,4)<=probe_res(nP,7),:);
        temp_behav_SW=[temp_test_res(:,1) nP*ones(size(temp_test_res,1),1) temp_test_res(:,4) temp_test_res(:,5)==3 temp_test_res(:,8)-probe_res(nP,3) temp_test_res(:,11)-temp_test_res(:,8) 1-temp_test_res(:,12:13) zeros(size(temp_test_res,1),2)];
        temp_behav_SW(temp_behav_SW(:,5)<-window_before_probes,:)=[];

        temp_SW=slow_Waves(slow_Waves(:,2)==nP,:);
        onset_Waves_FromProbes=temp_SW(:,5)/Fs+EEG.times(1)/1000;
        all_ElecFlags=nan(size(temp_behav_SW,1)-1,length(ChanLabels));
        for nT=1:size(temp_behav_SW,1)-1
            these_SW=temp_SW(onset_Waves_FromProbes>temp_behav_SW(nT,5) & onset_Waves_FromProbes<temp_behav_SW(nT+1,5),:);
            ElecFlags=hist(these_SW(:,3),1:length(ChanLabels));
            all_ElecFlags(nT,:)=ElecFlags>0;
        end
        temp_table_behav_SW=array2table([temp_behav_SW(1:size(temp_behav_SW,1)-1,:) all_ElecFlags],'VariableNames',[{'Block','Probe','Trial','Condition','DistToProbe','RT','FA','Miss','MS','Vigilance'},ChanLabels]);
        temp_table_behav_SW.MS=probe_res(nP,17)*ones(size(temp_table_behav_SW,1),1);
        temp_table_behav_SW.Vigilance=probe_res(nP,19)*ones(size(temp_table_behav_SW,1),1);
        table_behav_SW=[table_behav_SW ; temp_table_behav_SW];
    end
    table_behav_SW.SubID=repmat({SubID},size(table_behav_SW,1),1);
    all_table_behav_SW=[all_table_behav_SW ; table_behav_SW];
end
all_table_behav_SW.SubID=categorical(all_table_behav_SW.SubID);
all_table_behav_SW.MS=categorical(all_table_behav_SW.MS);
all_table_behav_SW.MS(all_table_behav_SW.MS=='1')='ON';
all_table_behav_SW.MS(all_table_behav_SW.MS=='2')='MW';
all_table_behav_SW.MS(all_table_behav_SW.MS=='3')='MB';
all_table_behav_SW.MS(all_table_behav_SW.MS=='4')='DK';
all_table_behav_SW.MS=removecats(all_table_behav_SW.MS);
all_table_behav_SW.Condition=categorical(all_table_behav_SW.Condition);
all_table_behav_SW.Condition(all_table_behav_SW.Condition=='0')='Go';
all_table_behav_SW.Condition(all_table_behav_SW.Condition=='1')='NoGo';
all_table_behav_SW.Condition=removecats(all_table_behav_SW.Condition);
writetable(all_table_behav_SW,[save_path filesep 'MW_MRI_SW_delta_Behav_PerTrial.txt']);

%%
mdl_FA=fitglme(all_table_behav_SW(all_table_behav_SW.Condition=='NoGo',:),'FA~1+Probe+MS+Vigilance+(1|SubID)','Distribution','Binomial');
mdl_Miss=fitglme(all_table_behav_SW(all_table_behav_SW.Condition=='Go' & (all_table_behav_SW.RT>0.15 | isnan(all_table_behav_SW.RT)),:),'Miss~1+Probe+MS+Vigilance+(1|SubID)','Distribution','Binomial');
mdl_RT=fitlme(all_table_behav_SW(all_table_behav_SW.Condition=='Go' & (all_table_behav_SW.RT>0.15 | isnan(all_table_behav_SW.RT)),:),'RT~1+Probe+MS+Vigilance+(1|SubID)');

%%
addpath((path_fieldtrip))
ft_defaults;

% ChanLabels={EEG.chanlocs.labels};
ChanLabels(find(ismember(ChanLabels,'FPz')))={'Fpz'};
ChanLabels(find(ismember(ChanLabels,'FP1')))={'Fp1'};
table_SW.Elec(find(ismember(table_SW.Elec,'FPz')))={'Fpz'};
table_SW.Elec(find(ismember(table_SW.Elec,'FP1')))={'Fp1'};
all_table_behav_SW = renamevars(all_table_behav_SW,["FPz","FP1"],["Fpz","Fp1"]);

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


%% Overall correlation between SW distibution and performance/mindstate  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FIG 5
% Correlation with behaviour
all_table_behav_SW(all_table_behav_SW.DistToProbe<-20,:)=[];
clear *_effect
newlabels=layout.label(1:end-2);
    fprintf('%2.0f/%2.0f\n',0,length(newlabels))
for nCh=1:length(newlabels)
    fprintf('\b\b\b\b\b\b%2.0f/%2.0f\n',nCh,length(newlabels))
    all_table_behav_SW.SW=all_table_behav_SW.(newlabels{nCh});
    mdl_FA_SW=fitglme(all_table_behav_SW(all_table_behav_SW.Condition=='NoGo' & (all_table_behav_SW.RT>0.15 | isnan(all_table_behav_SW.RT)),:),'FA~1+Block+SW+(1+Block|SubID)','Distribution','Binomial');
    mdl_Miss_SW=fitglme(all_table_behav_SW(all_table_behav_SW.Condition=='Go' & (all_table_behav_SW.RT>0.15 | isnan(all_table_behav_SW.RT)),:),'Miss~1+Block+SW+(1+Block|SubID)','Distribution','Binomial');
    mdl_RT_SW=fitlme(all_table_behav_SW(all_table_behav_SW.Condition=='Go'& (all_table_behav_SW.RT>0.15 | isnan(all_table_behav_SW.RT)),:),'RT~1+Block+SW+(1+Block|SubID)');


    %     mdl_FA=fitlme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_FA~1+Block+SW_density+(1+Block|SubID)');
    FA_effect(nCh,1)=mdl_FA_SW.Coefficients.tStat(match_str(mdl_FA_SW.CoefficientNames,'SW'));
    FA_effect(nCh,2)=mdl_FA_SW.Coefficients.pValue(match_str(mdl_FA_SW.CoefficientNames,'SW'));

    Miss_effect(nCh,1)=mdl_Miss_SW.Coefficients.tStat(match_str(mdl_Miss_SW.CoefficientNames,'SW'));
    Miss_effect(nCh,2)=mdl_Miss_SW.Coefficients.pValue(match_str(mdl_Miss_SW.CoefficientNames,'SW'));

    RT_effect(nCh,1)=mdl_RT_SW.Coefficients.tStat(match_str(mdl_RT_SW.CoefficientNames,'SW'));
    RT_effect(nCh,2)=mdl_RT_SW.Coefficients.pValue(match_str(mdl_RT_SW.CoefficientNames,'SW'));
end


%% Figure

% effect = [Miss_effect(:,1), FA_effect(:,1), RT_effect(:,1), MB_effect(:,1), MW_effect(:,1), VIG_effect(:,1)];
all_pV = [FA_effect(:,2); Miss_effect(:,2); RT_effect(:,2)];
FDR_Thr=0.05; %fdr(all_pV,0.05);
% cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
% minmax = ceil(max(max(abs(effect))));
% 
% %
figure;
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
subplot(1,3,1)
simpleTopoPlot_ft(FA_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(FA_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colormap(cmap2);
caxis([-1 1]*5)
title('FA', 'FontSize', 16)

subplot(1,3,2)
simpleTopoPlot_ft(Miss_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(Miss_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colormap(cmap2);
caxis([-1 1]*5)
title('Miss', 'FontSize', 16)

subplot(1,3,3)
simpleTopoPlot_ft(RT_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RT_effect(:,2)<0.05), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colormap(cmap2);
caxis([-1 1]*5)
title('Hit RT', 'FontSize', 16)

