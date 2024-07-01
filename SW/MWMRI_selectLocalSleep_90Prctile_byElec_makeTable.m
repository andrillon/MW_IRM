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
table_SW=array2table(zeros(0,15),'VariableNames',{'SubID','Block','Probe','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SART_Miss','SART_FA','SART_HitRT','SART_ON','SART_State','SART_Vig'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.Elec=categorical(table_SW.Elec);
table_SW.SART_State=categorical(table_SW.SART_State);

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
    if nFc==1
        addpath(genpath(path_eeglab));
    EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
    rmpath(genpath(path_eeglab));
    end
    
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
    
    for nP=unique(slow_Waves(:,2))'
        slow_Waves_perE=[];
        duration_of_probe=window_before_probes/60;
        
        temp_test_res=test_res(test_res(:,1)==probe_res(nP,5) & test_res(:,4)<=probe_res(nP,7),:);
        temp_go=temp_test_res(temp_test_res(:,5)~=3,:);
             temp_nogo=temp_test_res(temp_test_res(:,5)==3,:);
   if size(temp_go,1)<18 || size(temp_nogo,1)<2
            warning(sprintf('... skipping probe %g in block %g from %s because missing data',probe_res(nP,1),probe_res(nP,5),SubID))
            continue;
        end
        temp_go=temp_go(end-17:end,:);
        temp_nogo=temp_nogo(end-1:end,:);
        
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
        
        table_SW.SART_Miss(table_length+(1:length(ChanLabels)))=repmat(Miss,length(ChanLabels),1);
        table_SW.SART_HitRT(table_length+(1:length(ChanLabels)))=repmat(HitRT,length(ChanLabels),1);
        table_SW.SART_FA(table_length+(1:length(ChanLabels)))=repmat(FA,length(ChanLabels),1);
        table_SW.SART_State(table_length+(1:length(ChanLabels)))=repmat({this_State},length(ChanLabels),1);
        table_SW.SART_Vig(table_length+(1:length(ChanLabels)))=repmat(this_Vig,length(ChanLabels),1);
        
    end
end
writetable(table_SW,[save_path filesep 'MW_MRI_SW_Behav_PerProbe.txt']);
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


%% Overall correlation between SW distibution and performance/mindstate  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FIG 5
% Correlation with behaviour
clear *_effect
newlabels=layout.label(1:end-2);
for nCh=1:length(newlabels)
        fprintf('channel %s %g/%g\n',newlabels{nCh},nCh,length(newlabels))

    mdl_FA=fitlme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_FA~1+Block+SW_density+(1|SubID)');
    FA_effect(nCh,1)=mdl_FA.Coefficients.tStat(match_str(mdl_FA.CoefficientNames,'SW_density'));
    FA_effect(nCh,2)=mdl_FA.Coefficients.pValue(match_str(mdl_FA.CoefficientNames,'SW_density'));
    
    mdl_Miss=fitlme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_Miss~1+Block+SW_density+(1|SubID)');
    Miss_effect(nCh,1)=mdl_Miss.Coefficients.tStat(match_str(mdl_Miss.CoefficientNames,'SW_density'));
    Miss_effect(nCh,2)=mdl_Miss.Coefficients.pValue(match_str(mdl_Miss.CoefficientNames,'SW_density'));
    
    mdl_RT=fitlme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_HitRT~1+Block+SW_density+(1|SubID)');
    RT_effect(nCh,1)=mdl_RT.Coefficients.tStat(match_str(mdl_RT.CoefficientNames,'SW_density'));
    RT_effect(nCh,2)=mdl_RT.Coefficients.pValue(match_str(mdl_RT.CoefficientNames,'SW_density'));
    
end



% Correlation with Mindstate
%%
%table2_SW.SART_MB2=table2_SW.SART_MB+table2_SW.SART_DK;
% table_SW.SART_MW(table_SW.SART_MW~=1 & table_SW.SART_ON~=1)=NaN;
% table_SW.SART_MB(table_SW.SART_MB~=1 & table_SW.SART_ON~=1)=NaN;
table_SW.SART_ON=nan(size(table_SW,1),1);
table_SW.SART_ON(table_SW.SART_State=='ON')=1;
table_SW.SART_ON(table_SW.SART_State~='ON')=0;

table_SW.SART_MW=nan(size(table_SW,1),1);
table_SW.SART_MW(table_SW.SART_State=='MW')=1;
table_SW.SART_MW(table_SW.SART_State=='ON')=0;

table_SW.SART_MB=nan(size(table_SW,1),1);
table_SW.SART_MB(table_SW.SART_State=='MB' | table_SW.SART_State=='DK')=1;
table_SW.SART_MB(table_SW.SART_State=='ON')=0;

table_SW.SART_Vig2=nan(size(table_SW,1),1);
table_SW.SART_Vig2(table_SW.SART_Vig==1 | table_SW.SART_Vig==2)=0;
table_SW.SART_Vig2(table_SW.SART_Vig==3 | table_SW.SART_Vig==4)=1;

for nCh=1:length(newlabels)
    fprintf('channel %s %g/%g\n',newlabels{nCh},nCh,length(newlabels))
    mdl_MW=fitglme(table_SW(table_SW.Elec==newlabels{nCh} & isnan(table_SW.SART_MW)==0,:),'SART_MW~1+Block+SW_density+(1+Block|SubID)','Distribution','binomial');
    MW_effect(nCh,1)=mdl_MW.Coefficients.tStat(match_str(mdl_MW.CoefficientNames,'SW_density'));
    MW_effect(nCh,2)=mdl_MW.Coefficients.pValue(match_str(mdl_MW.CoefficientNames,'SW_density'));
    
    mdl_MB=fitglme(table_SW(table_SW.Elec==newlabels{nCh} & isnan(table_SW.SART_MB)==0,:),'SART_MB~1+Block+SW_density+(1+Block|SubID)','Distribution','binomial');
    MB_effect(nCh,1)=mdl_MB.Coefficients.tStat(match_str(mdl_MB.CoefficientNames,'SW_density'));
    MB_effect(nCh,2)=mdl_MB.Coefficients.pValue(match_str(mdl_MB.CoefficientNames,'SW_density'));
    
    mdl_ON=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_ON~1+Block+SW_density+(1+Block|SubID)','Distribution','binomial');
    ON_effect(nCh,1)=mdl_ON.Coefficients.tStat(match_str(mdl_ON.CoefficientNames,'SW_density'));
    ON_effect(nCh,2)=mdl_ON.Coefficients.pValue(match_str(mdl_ON.CoefficientNames,'SW_density'));
    
    mdl_VIG=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'SART_Vig2~1+Block+SW_density+(1+Block|SubID)');
    VIG_effect(nCh,1)=mdl_VIG.Coefficients.tStat(match_str(mdl_VIG.CoefficientNames,'SW_density'));
    VIG_effect(nCh,2)=mdl_VIG.Coefficients.pValue(match_str(mdl_VIG.CoefficientNames,'SW_density'));

end

%% Figure

effect = [Miss_effect(:,1), FA_effect(:,1), RT_effect(:,1), MB_effect(:,1), MW_effect(:,1), VIG_effect(:,1)];
all_pV = [Miss_effect(:,2); FA_effect(:,2); RT_effect(:,2); MB_effect(:,2); MW_effect(:,2); VIG_effect(:,2)];
FDR_Thr=0.05; %fdr(all_pV,0.05);
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
minmax = ceil(max(max(abs(effect))));

%
figure;

subplot(2,3,1)
simpleTopoPlot_ft(Miss_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(Miss_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('Misses', 'FontSize', 16)

subplot(2,3,2);
simpleTopoPlot_ft(FA_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(FA_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('False Alarms', 'FontSize', 16)

subplot(2,3,4);
simpleTopoPlot_ft(-ON_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON', 'FontSize', 16)

subplot(2,3,3);
simpleTopoPlot_ft(VIG_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(VIG_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('Vigilance', 'FontSize', 16)

subplot(2,3,5);
simpleTopoPlot_ft(MW_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(MW_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('Mind Wandering', 'FontSize', 16)

subplot(2,3,6);
simpleTopoPlot_ft(MB_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(MB_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('Mind Blanking', 'FontSize', 16)
% c = colorbar; c.Label.String = 't-value'; c.Label.FontSize = 14; c.Label.Rotation = 270; c.Label.Position(1) = 2.5; c.Ticks = [-8 8]; c.FontSize = 14;

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

