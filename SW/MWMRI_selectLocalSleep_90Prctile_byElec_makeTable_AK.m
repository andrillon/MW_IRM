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
table_SW=array2table(zeros(0,19),'VariableNames',{'SubID','Block','Probe','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SART_Miss','SART_FA','SART_HitRT','SART_ON','SART_State','SART_Vig', 'ROI1', 'ROI2', 'ROI3', 'ROI4'});
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
    ROI1 = importdata([ROIpath 'Preprobe10sARV_not_ot03422_s' SubID(end-2:end) '.mat']); % Bilat Ant cing
    ROI2 = importdata([ROIpath 'Preprobe10sARV_not_ot2-3830_s' SubID(end-2:end) '.mat']); % Bilat mid-cing, L post cing
    ROI3 = importdata([ROIpath 'Preprobe10sARV_not_ot_14-6430_s' SubID(end-2:end) '.mat']); %R cuenous, precuneous
    ROI4 = importdata([ROIpath 'Preprobe10sARV_ot_not_36-42-4_s' SubID(end-2:end) '.mat']); % FFA

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

        dataProbeind = unique(slow_Waves(slow_Waves(:,2)==nP,16));

        this_ROI1=ROI1(dataProbeind,:);
        table_SW.ROI1(table_length+(1:length(ChanLabels)))=repmat(this_ROI1,length(ChanLabels),1);
        this_ROI2=ROI2(dataProbeind,:);
        table_SW.ROI2(table_length+(1:length(ChanLabels)))=repmat(this_ROI2,length(ChanLabels),1);
        this_ROI3=ROI3(dataProbeind,:);
        table_SW.ROI3(table_length+(1:length(ChanLabels)))=repmat(this_ROI3,length(ChanLabels),1);
        this_ROI4=ROI4(dataProbeind,:);
        table_SW.ROI4(table_length+(1:length(ChanLabels)))=repmat(this_ROI4,length(ChanLabels),1);

        
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


%% Figure

effect = [ON_D_effect(:,1), ON_A_effect(:,1), ON_F_effect(:,1), ON_Do_effect(:,1), ON_U_effect(:,1)];
%all_pV = [ON_D_effect(:,2); ON_A_effect(:,2); ON_F_effect(:,2); ON_Do_effect(:,2); ON_U_effect(:,2)];
all_pV = [ON_D_effect(:,2)];

% effect = [ROI1_effect(:,1), ROI1_effect(:,3), ROI1_effect(:,5)];
% all_pV = [ROI1_effect(:,2), ROI1_effect(:,4), ROI1_effect(:,6)];
FDR_Thr=0.05;%fdr(all_pV,0.05);
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
minmax = ceil(max(max(abs(effect))));

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



%%



for nCh=1:length(newlabels)
    mdl_ROI1_D=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROI2~1+Block+SW_density+(1+Block|SubID)','Distribution','Normal');
    ROI_effect(nCh,1)=mdl_ROI1_D.Coefficients.tStat(match_str(mdl_ROI1_D.CoefficientNames,'SW_density'));
    ROI_effect(nCh,2)=mdl_ROI1_D.Coefficients.pValue(match_str(mdl_ROI1_D.CoefficientNames,'SW_density'));
    fprintf('Electrode %s\n',newlabels{nCh})

    mdl_ROI1_A=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROI2~1+Block+SW_amplitude+(1+Block|SubID)','Distribution','Normal');
    ROI_A_effect(nCh,1)=mdl_ROI1_A.Coefficients.tStat(match_str(mdl_ROI1_A.CoefficientNames,'SW_amplitude'));
    ROI_A_effect(nCh,2)=mdl_ROI1_A.Coefficients.pValue(match_str(mdl_ROI1_A.CoefficientNames,'SW_amplitude'));

    mdl_ROI1_F=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROI2~1+Block+SW_frequency+(1+Block|SubID)','Distribution','Normal');
    ROI_F_effect(nCh,1)=mdl_ROI1_F.Coefficients.tStat(match_str(mdl_ROI1_F.CoefficientNames,'SW_frequency'));
    ROI_F_effect(nCh,2)=mdl_ROI1_F.Coefficients.pValue(match_str(mdl_ROI1_F.CoefficientNames,'SW_frequency'));

    mdl_ROI1_Do=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROI2~1+Block+SW_downslope+(1+Block|SubID)','Distribution','Normal');
    ROI_Do_effect(nCh,1)=mdl_ROI1_Do.Coefficients.tStat(match_str(mdl_ROI1_Do.CoefficientNames,'SW_downslope'));
    ROI_Do_effect(nCh,2)=mdl_ROI1_Do.Coefficients.pValue(match_str(mdl_ROI1_Do.CoefficientNames,'SW_downslope'));

    mdl_ROI1_U=fitglme(table_SW(table_SW.Elec==newlabels{nCh},:),'ROI2~1+Block+SW_upslope+(1+Block|SubID)','Distribution','Normal');
    ROI_U_effect(nCh,1)=mdl_ROI1_U.Coefficients.tStat(match_str(mdl_ROI1_U.CoefficientNames,'SW_upslope'));
    ROI_U_effect(nCh,2)=mdl_ROI1_U.Coefficients.pValue(match_str(mdl_ROI1_U.CoefficientNames,'SW_upslope'));

end


effect = [ROI_effect(:,1), ROI_A_effect(:,1), ROI_F_effect(:,1), ROI_Do_effect(:,1), ROI_U_effect(:,1) ];
FDR_Thr=0.05; %fdr(all_pV,0.05);
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);
minmax = ceil(max(max(abs(effect))));


figure;
set(gcf, 'Position', [100, 100, 800, 600]); % Set the figure position and size


subplot(2,2,1)
simpleTopoPlot_ft(ROI_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ROI_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI-1 (right Fusform.) - Density', 'FontSize', 16)

subplot(2,2,2)
simpleTopoPlot_ft(ROI_A_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ROI_A_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI-1 (right Fusform.) - Amplitude', 'FontSize', 16)

subplot(2,2,3)
simpleTopoPlot_ft(ROI_Do_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ROI_Do_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI-1 (right Fusform.) - Downslope', 'FontSize', 16)

subplot(2,2,4)
simpleTopoPlot_ft(ROI_U_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ROI_U_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('ROI-1 (right Fusform.) - Upslope', 'FontSize', 16)



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

