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
    
    Fs=500;
    
    load([save_path filesep 'prct_10s_DSS_SW_' SubID]);%,'slow_Waves','paramSW','ChanLabels')
    
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
%%

% select relevant files, here baseline blocks
files=dir([data_path filesep filesep 'RestingState/*clean_nt.set']);
table_SWR=array2table(zeros(0,9),'VariableNames',{'SubID','Block','Probe','Elec','SW_density','SW_amplitude','SW_frequency','SW_downslope','SW_upslope'});
table_SWR.SubID=categorical(table_SWR.SubID);
table_SWR.Elec=categorical(table_SWR.Elec);

%% loop across trials for baseline blocks
nFc=0;
all_SW_probes=[];
window_before_probes=30; % in seconds
for nF=1:length(files)
    % load file with EEGlab
    fprintf('... file: %s\n',files(nF).name)
    
    SubID=files(nF).name;
    sep=findstr(SubID,'_RS_clean_nt.set'); %
    
    if isempty(sep)
        SubID=SubID(1:end-9);
    else
        SubID=SubID(1:sep(1)-1);
    end
    if exist([save_path filesep 'DSS_RS_allSW_' SubID '.mat'])==0  %ICA_RS_allSW_
        continue;
    end
   
    % load EEG
    addpath(genpath(path_eeglab));
    EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);

    rmpath(genpath(path_eeglab));

    ChanLabels={EEG.chanlocs.labels};
    EEG.data=EEG.data(~ismember(ChanLabels,{'ECG','ECG2'}),:,:);
    ChanLabels=ChanLabels(~ismember(ChanLabels,{'ECG','ECG2'}));
    nFc=nFc+1;

    load([save_path filesep 'prct_DSS_RS_SW_' SubID]);%,'slow_Waves','paramSW','ChanLabels')

     for nP=unique(slow_Waves(:,2))'
        slow_Waves_perE=[];
        duration_of_probe=window_before_probes/60;
        

            for nE=1:length(ChanLabels)
                    slow_Waves_perE=[slow_Waves_perE ; [sum(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP)/duration_of_probe nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,4)) nanmean(1./((slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,7)-slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,5))/Fs)) ...
                        nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,12)) nanmean(slow_Waves(slow_Waves(:,3)==nE & slow_Waves(:,2)==nP,13)) NaN]];
            end
            table_length=size(table_SWR,1);
            
            table_SWR.SubID(table_length+(1:length(ChanLabels)))=repmat({SubID},length(ChanLabels),1);
            table_SWR.Block(table_length+(1:length(ChanLabels)))=repmat(1,length(ChanLabels),1);
            table_SWR.Probe(table_length+(1:length(ChanLabels)))=repmat(nP,length(ChanLabels),1);
            table_SWR.Elec(table_length+(1:length(ChanLabels)))=ChanLabels;
            
            table_SWR.SW_density(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,1);
            table_SWR.SW_amplitude(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,2);
            table_SWR.SW_frequency(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,3);
            table_SWR.SW_downslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,4);
            table_SWR.SW_upslope(table_length+(1:length(ChanLabels)))=slow_Waves_perE(:,5);
      end
    
end


%%

% Filter for Block 1
table_SWR_Block1 = table_SWR(table_SWR.Block == 1, :);
table_SW_Block1 = table_SW(:, :);


% Average data across probes for each SubID and Elec
varsToAverage = {'SW_density', 'SW_amplitude', 'SW_frequency', 'SW_downslope', 'SW_upslope'};
groupVars = {'SubID', 'Elec'};
table_SW_Block1_avg = varfun(@nanmean, table_SW_Block1, 'InputVariables', varsToAverage, 'GroupingVariables', groupVars);
table_SWR_Block1_avg = varfun(@nanmean, table_SWR_Block1, 'InputVariables', varsToAverage, 'GroupingVariables', groupVars);


% Merge the average tables by SubID and Elec
mergedTable = outerjoin(table_SW_Block1_avg, table_SWR_Block1_avg, 'Keys', {'SubID', 'Elec'}, 'MergeKeys', true);
% Define new names
newNames = {'SubID', 'Elec', 'GroupCount_SW', 'SW_Density_SW', 'SW_Amplitude_SW', ...
            'SW_Frequency_SW', 'SW_Downslope_SW', 'SW_Upslope_SW', 'GroupCount_SWR', ...
            'SW_Density_SWR', 'SW_Amplitude_SWR', 'SW_Frequency_SWR', 'SW_Downslope_SWR', 'SW_Upslope_SWR'};

% Rename variables in the merged table
mergedTable = renamevars(mergedTable, ...
    {'SubID', 'Elec', 'GroupCount_table_SW_Block1_avg', 'nanmean_SW_density_table_SW_Block1_avg', ...
     'nanmean_SW_amplitude_table_SW_Block1_avg', 'nanmean_SW_frequency_table_SW_Block1_avg', ...
     'nanmean_SW_downslope_table_SW_Block1_avg', 'nanmean_SW_upslope_table_SW_Block1_avg', ...
     'GroupCount_table_SWR_Block1_avg', 'nanmean_SW_density_table_SWR_Block1_avg', ...
     'nanmean_SW_amplitude_table_SWR_Block1_avg', 'nanmean_SW_frequency_table_SWR_Block1_avg', ...
     'nanmean_SW_downslope_table_SWR_Block1_avg', 'nanmean_SW_upslope_table_SWR_Block1_avg'}, ...
    newNames);


% Find common SubIDs
uniqueSubIDs = intersect(unique(table_SW_Block1_avg.SubID), unique(table_SWR_Block1_avg.SubID));

% Loop through each SubID to correlate and plot
for i = 1:length(uniqueSubIDs)
    subID = string(uniqueSubIDs(i));
   
    % Filter merged data for the current SubID
    subData = mergedTable(strcmp(string(mergedTable.SubID), subID), :);

    % Calculate correlation for nanmean_SW_density
    [r, p] = corr(subData.SW_Density_SW, subData.SW_Density_SWR,'Type','Pearson','Rows','complete'); 
    % Calculate coefficients for the line of best fit
    coeffs = polyfit(subData.SW_Density_SW, subData.SW_Density_SWR, 1);  % Linear fit (degree 1)
    slope = coeffs(1);
    intercept = coeffs(2);


    % Create scatter plot for the current SubID
    figure;
    scatter(subData.SW_Density_SW, subData.SW_Density_SWR, 'filled');
    xlim([0,10]); 
    ylim([0,10]); 
    hold on;
    % Plot line of best fit
    xFit = xlim();  % Get the current x-axis limits for the fit line
    yFit = polyval(coeffs, xFit);  % Generate y data based on the fit
    plot(xFit, yFit, 'r-', 'LineWidth', 1);  % Plot the line of best fit in red
    title(['SW Density ', subID,sprintf(' (r=%.2f, p=%.3f)', r, p)]);
    xlabel('SW Task (All blocks)');
    ylabel('SW Resting');
    grid on;
    format_fig;
    hold off;  

    % Save the figure
    outputDir = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/SW_RestVSTask/';
    filename = sprintf('%s%s_RSvsTaskBlockall.png', outputDir, subID);  % Construct filename
    saveas(gcf, filename);  % Save the figure as a PNG file

    % Close the figure to free up system resources and avoid display overload
    close(gcf);

end


% Initialize the correlation table
corrTable = table();

% Loop over all variable pairs and compute correlation
for i = 1:length(uniqueSubIDs)
    subID = string(uniqueSubIDs(i));
     % Filter merged data for the current SubID
    subData = mergedTable(strcmp(string(mergedTable.SubID), subID), :);

    % Calculate correlation for nanmean_SW_density
    [r, p] = corr(subData.SW_Density_SW, subData.SW_Density_SWR,'Type','Pearson','Rows','complete'); 

    % Store results in the table
    corrTable.subID(i) = {subID};  % Store r value
    corrTable.Coef__Denisty(i) = {r};  % Store r value
    corrTable.Pval__Denisty(i) = {p}; % Store p value

    [r, p] = corr(subData.SW_Amplitude_SW, subData.SW_Amplitude_SWR,'Type','Pearson','Rows','complete'); 
    corrTable.Coef__Amp(i) = {r};  % Store r value
    corrTable.Pval__Amp(i) = {p}; % Store p value

    [r, p] = corr(subData.SW_Frequency_SW, subData.SW_Frequency_SWR,'Type','Pearson','Rows','complete'); 
    corrTable.Coef__Freq(i) = {r};  % Store r value
    corrTable.Pval__Freq(i) = {p}; % Store p value

    [r, p] = corr(subData.SW_Downslope_SW, subData.SW_Downslope_SWR,'Type','Pearson','Rows','complete'); 
    corrTable.Coef__Dslope(i) = {r};  % Store r value
    corrTable.Pval__Dslope(i) = {p}; % Store p value

    [r, p] = corr(subData.SW_Upslope_SW, subData.SW_Upslope_SWR,'Type','Pearson','Rows','complete'); 
    corrTable.Coef__Uslope(i) = {r};  % Store r value
    corrTable.Pval__Uslope(i) = {p}; % Store p value

end




for i = 2:11 % Start from 2 to skip 'SubID'
    figure; % Open a new figure for each histogram
    datatoplot = cell2mat(corrTable{:,i});
    histogram(datatoplot, 'NumBins', 15, 'FaceColor', 'blue'); % Plot histogram of the i-th column
    title(['Histogram of ', corrTable.Properties.VariableNames{i}]); % Set title dynamically based on column name
    xlabel(['Values of ', corrTable.Properties.VariableNames{i}]); % Label x-axis
    ylabel('Frequency'); % Label y-axis
    grid on;
    format_fig;
end

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

figure;
simpleTopoPlot_ft(ON_effect(:,1), layout,'on',[],0,1);
ft_plot_lay_me(layout, 'chanindx', find(ON_effect(:,2)<FDR_Thr), 'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','yes')
%colorbar;
colormap(cmap2);
caxis([-1 1]*5)
title('OFF vs ON', 'FontSize', 16)
c = colorbar; c.Label.String = 't-value'; c.Label.FontSize = 14; c.Label.Rotation = 270; c.Label.Position(1) = 4; c.FontSize = 14;



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

