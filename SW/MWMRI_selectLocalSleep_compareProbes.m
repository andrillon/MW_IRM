%%%%% Preprocessing check
%%%%% Basic analyses on the preprocessed, filtered and epoch EEG data

%% Init
clear all;
close all;

run ../localdef.m

% adding relevant toolboxes to the path
addpath(genpath(lscpTools_path))
addpath((path_fieldtrip))
ft_defaults
addpath(genpath(path_eeglab));

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'MWMRI*clean.set']);

myERP_Elec={'Fz','Cz','Pz','Oz'};
States={'ON','MW','MB','??'};
%%
SW_table=array2table(zeros(0,15),'VariableNames',{'SubID','Block','Probe','Elec','State','Vigilance','Misses','FAs','Hit_RT','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope','SW_threshold'});
SW_table.SubID=categorical(SW_table.SubID);
SW_table.Elec=categorical(SW_table.Elec);
SW_table.State=categorical(SW_table.State);


%% loop across trials for baseline blocks
redo=1;
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
    if exist([save_path filesep 'prct_ICA_SW_' SubID '.mat'])==0
        continue;
    end

    % load behaviour
    file_behav=dir([data_path filesep '..' filesep '..' filesep 'Behav' filesep 'wanderIM_behavres_s' SubID(6:end) '*.mat']);
    if ~isempty(file_behav)
        load([file_behav.folder filesep file_behav.name])
    else
        continue;
    end
    

    %      size(EEG.data)
    nFc=nFc+1;

    load([save_path filesep 'prct_ICA_SW_' SubID])
    %     slow_Waves(slow_Waves(:,3)==length(ChanLabels),:)=[];
    
    for nBl=1:max(slow_Waves(:,2))
        slow_Waves_perE=[];
        duration_probe=25/60;
        Fs=500;
        for nE=1:length(ChanLabels)
            slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl)/duration_probe nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,5))/Fs)) ...
                nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,13)) NaN nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,9)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nBl,11))]];
        end
        
        table_length=size(SW_table,1);
        SW_table.SubID(table_length+(1:length(ChanLabels)))=repmat({SubID},length(ChanLabels),1);
        SW_table.Elec(table_length+(1:length(ChanLabels)))=ChanLabels;
        SW_table.State(table_length+(1:length(ChanLabels)))=repmat(States(probe_res(nBl,17)),length(ChanLabels),1);
        
        SW_table.Probe(table_length+(1:length(ChanLabels)))=repmat(nBl,length(ChanLabels),1);
        SW_table.Vigilance(table_length+(1:length(ChanLabels)))=repmat(probe_res(nBl,19),length(ChanLabels),1);
        SW_table.Block(table_length+(1:length(ChanLabels)))=repmat(probe_res(nBl,5),length(ChanLabels),1);
        
        SW_table.SW_density(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,1);
        SW_table.SW_amplitude(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,2);
        SW_table.SW_frequency(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,3);
        SW_table.SW_downslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,4);
        SW_table.SW_upslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,5);
        SW_table.SW_threshold(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,6);
        
        temp_test_res=test_res(test_res(:,1)==probe_res(nBl,5) & test_res(:,4)>=probe_res(nBl,7)-25 & test_res(:,4)<probe_res(nBl,7),:);
        FA=nanmean(temp_test_res(:,12)==0);
        Misses=nanmean(temp_test_res(:,13)==0);
        Hit_RT=nanmean(temp_test_res(temp_test_res(:,5)~=3,11)-temp_test_res(temp_test_res(:,5)~=3,8));
        SW_table.Misses(table_length+(1:length(ChanLabels)))=repmat(Misses,length(ChanLabels),1);
        SW_table.FAs(table_length+(1:length(ChanLabels)))=repmat(FA,length(ChanLabels),1);
        SW_table.Hit_RT(table_length+(1:length(ChanLabels)))=repmat(Hit_RT,length(ChanLabels),1);
    end
end


%%
figure;
temp_plot=[];
for nP=1:40
    temp_plot(nP)=mean(SW_table.SW_density(SW_table.Elec=='Cz' & SW_table.Probe==nP));
end
plot(1:40,temp_plot);
hold on; format_fig;
for k=[10.5 20.5 30.5]
    line([1 1]*k,ylim,'Color','k','LineStyle','--') 
end


%% Topography
SW_table.Elec(SW_table.Elec=='FPz')='Fpz';
SW_table.Elec(SW_table.Elec=='FP1')='Fp1';


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

%%
figure
topo_plot=[];
for nCh=1:length(layout.label)-2
    topo_plot(nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh})); %squeeze(mean(mean(SW_dens_perProbe,2),1));
end
simpleTopoPlot_ft(topo_plot', layout,'labels',[],0,1);

%%
caxis_all=[];
figure;
for k=1:4
    subplot(2,2,k)
    topo_plot=[];
    for nCh=1:length(layout.label)-2
        topo_plot(nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.Vigilance==k)); %squeeze(mean(mean(SW_dens_perProbe,2),1));
    end
    simpleTopoPlot_ft(topo_plot', layout,'on',[],0,1);
    colorbar;
        title(sprintf('VIG: %g',k))
        caxis_all=[caxis_all ; [min(topo_plot) max(topo_plot)]];

end
for k=1:4
    subplot(2,2,k)
    caxis([min(caxis_all(:,1)) max(caxis_all(:,2))])
end
%%
figure;
titleNames={'ON','MW','MB','??'};
caxis_all=[];
for k=1:4
    subplot(2,2,k)
    topo_plot=[];
    for nCh=1:length(layout.label)-2
        topo_plot(nCh)=nanmean(SW_table.SW_density(SW_table.Elec==layout.label{nCh} & SW_table.State==States{k})); %squeeze(mean(mean(SW_dens_perProbe,2),1));
    end
    simpleTopoPlot_ft(topo_plot', layout,'on',[],0,1);
    colorbar;
    title(titleNames{k})
        caxis_all=[caxis_all ; [min(topo_plot) max(topo_plot)]];
end
for k=1:4
    subplot(2,2,k)
    caxis([min(caxis_all(:,1)) max(caxis_all(:,2))])
end
%%
figure;
topo_plot=[];
for nCh=1:length(layout.label)-2
    topo_plot(nCh)=nanmean(SW_table.SW_downslope(SW_table.Elec==layout.label{nCh})); %squeeze(mean(mean(SW_dens_perProbe,2),1));
end
simpleTopoPlot_ft(topo_plot', layout,'on',[],0,1);
colorbar;

%%
SW_table.Elec=categorical(SW_table.Elec);

for nCh=1:length(layout.label)-2
    
    mdl_FA=fitlm(SW_table(SW_table.Elec==layout.label{nCh},:),'FAs~1+Probe+SW_density');
    FA_effect(nCh,1)=mdl_FA.Coefficients.tStat(end);
    FA_effect(nCh,2)=mdl_FA.Coefficients.pValue(end);
    
    mdl_Miss=fitlm(SW_table(SW_table.Elec==layout.label{nCh},:),'Misses~1+Probe+SW_density');
    Miss_effect(nCh,1)=mdl_Miss.Coefficients.tStat(end);
    Miss_effect(nCh,2)=mdl_Miss.Coefficients.pValue(end);
    
    mdl_RT=fitlm(SW_table(SW_table.Elec==layout.label{nCh},:),'Hit_RT~1+Probe+SW_density');
    RT_effect(nCh,1)=mdl_RT.Coefficients.tStat(end);
    RT_effect(nCh,2)=mdl_RT.Coefficients.pValue(end);
    
end


%%
cmap2=cbrewer('div','RdBu',64); cmap2=flipud(cmap2);

figure;
subplot(1,3,2);
simpleTopoPlot_ft(Miss_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(Miss_effect(:,2)<fdr(Miss_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar; colormap(cmap2);
caxis([-1 1]*max(abs(Miss_effect(:,1))))
title('SW * Miss')

subplot(1,3,1);
simpleTopoPlot_ft(FA_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(FA_effect(:,2)<fdr(FA_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar; colormap(cmap2);
caxis([-1 1]*max(abs(FA_effect(:,1))))
title('SW * FA')

subplot(1,3,3);
simpleTopoPlot_ft(RT_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(RT_effect(:,2)<fdr(RT_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar; colormap(cmap2);
caxis([-1 1]*max(abs(RT_effect(:,1))))
title('SW * RT')

%%
SW_table.State2=SW_table.State;
SW_table.State2(SW_table.State2=='??')='MB';
SW_table.State2=removecats(SW_table.State2);

for nCh=1:length(layout.label)-2
    
    mdl_State=fitlm(SW_table(SW_table.Elec==layout.label{nCh},:),'SW_density~1+Probe+State2');
    MW_effect(nCh,1)=mdl_State.Coefficients.tStat(end);
    MW_effect(nCh,2)=mdl_State.Coefficients.pValue(end);
    
    MB_effect(nCh,1)=mdl_State.Coefficients.tStat(end-1);
    MB_effect(nCh,2)=mdl_State.Coefficients.pValue(end-1);
    
    mdl_Vigilance=fitlm(SW_table(SW_table.Elec==layout.label{nCh},:),'SW_density~1+Probe+Vigilance');
    Vig_effect(nCh,1)=mdl_Vigilance.Coefficients.tStat(end-1);
    Vig_effect(nCh,2)=mdl_Vigilance.Coefficients.pValue(end-1);
    
    
end

%%
figure;
subplot(1,3,1);
simpleTopoPlot_ft(MW_effect(:,1), layout,'on',[],0,1);
% ft_plot_lay_me(layout, 'chanindx', find(MW_effect(:,2)<fdr(MW_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar; colormap(cmap2);
caxis([-1 1]*max(abs(MW_effect(:,1))))
title('SW * MW')

subplot(1,3,2);
simpleTopoPlot_ft(MB_effect(:,1), layout,'on',[],0,1);
% ft_plot_lay_me(layout, 'chanindx', find(MB_effect(:,2)<fdr(MB_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar; colormap(cmap2);
caxis([-1 1]*max(abs(MB_effect(:,1))))
title('SW * MB')

subplot(1,3,3);
simpleTopoPlot_ft(Vig_effect(:,1), layout,'on',[],0,1);
% ft_plot_lay_me(layout, 'chanindx', find(Vig_effect(:,2)<fdr(Vig_effect(:,2),0.05)), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
colorbar; colormap(cmap2);
caxis([-1 1]*max(abs(Vig_effect(:,1))))
title('SW * Vig')
