 
%% Init

% clear all;
% close all;
cd /Users/kuszti/Documents/Git/MW_IRM/SW/
run ../localdef.m

files=dir([save_path filesep 'prct_DSS_SW_M*']);

for nF=1:length(files)

saveWaves =load([files(nF).folder filesep files(nF).name]);
slow_Waves= saveWaves.slow_Waves;


%%
% Change onset and end of wave time into seconds from the probe
Trialtimes = -25000:2:-2;
slow_Waves(:,5)=Trialtimes(slow_Waves(:,5));
slow_Waves(:,6)=Trialtimes(slow_Waves(:,6));
slow_Waves(:,7)=Trialtimes(slow_Waves(:,7));
slow_Waves(:,8)=Trialtimes(slow_Waves(:,8));
slow_Waves(:,10)=Trialtimes(slow_Waves(:,10));


% Select Slow waves based on time-window
window_before_probes=20; % in seconds
slow_Waves(slow_Waves(:,5)<-window_before_probes*1000,:)=[];

% box_waves=[];
% for nP=1:40
%     for nEl=unique(slow_Waves(:,3))'
%         temp_waves=slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,3)==nEl,:);
%         temp_vec=zeros(1,10*500);
%         for nW=1:size(temp_waves,1)
%             temp_vec((temp_waves(nW,5):temp_waves(nW,8))+10*500)=abs(temp_waves(nW,9));
%         end
%         box_waves(nP,nEl,:)=temp_vec;
%     end
% end
% 
% figure; plot (squeeze(box_waves(1,1,:)))
% hold on
% for e = 1:size(box_waves,2); plot (squeeze(box_waves(1,e,:))); end
% 
% figure; imagesc (mean(box_waves, 3)')

%%
% load fMRI timing
SubID=files(nF).name;
SubID=SubID(18:end-4);
subFolder = ['s' (SubID)];
fprintf('Running subject: %s \n',subFolder);


BlockDatasetP = dir([TimingPath subFolder '/*/*ProbeBlock*.csv']); 
BlockDataset = dir([TimingPath subFolder '/*/b*_Block*.csv']); 
Block = [];
BlockP = [];
for nbl = 1:size(BlockDataset,1)
    Block{nbl} = readmatrix([BlockDataset(nbl).folder filesep BlockDataset(nbl).name]);
    BlockP{nbl} = readmatrix([BlockDatasetP(nbl).folder filesep BlockDatasetP(nbl).name]);
end

%%
% convert into onset vector, duration vector, and parameric modulation
% based on amp.

% separet into blocks
Block1sws = slow_Waves(slow_Waves(:,2)>=1 & slow_Waves(:,2)<=10,:);
Block2sws = slow_Waves(slow_Waves(:,2)>=11 & slow_Waves(:,2)<=20,:);
Block3sws = slow_Waves(slow_Waves(:,2)>=21 & slow_Waves(:,2)<=30,:);
Block4sws = slow_Waves(slow_Waves(:,2)>=31 & slow_Waves(:,2)<=40,:);


% for s218
if strcmp(subFolder, 's218')==1
    Block1sws = slow_Waves(slow_Waves(:,2)>=1 & slow_Waves(:,2)<=17,:);
    Block2sws = slow_Waves(slow_Waves(:,2)>=19 & slow_Waves(:,2)<=20,:);
    tempBlock1 = BlockP{1};
    tempBlock1(11:17,1)=tempBlock1(11:17,1)+10;
    BlockP{1} = tempBlock1;
end

if strcmp(subFolder, 's228')==1
    Block3sws = slow_Waves(slow_Waves(:,2)>=21 & slow_Waves(:,2)<=38,:);
    Block4sws = slow_Waves(slow_Waves(:,2)==40,:);
    tempBlock3 = BlockP{3};
    tempBlock3(11:18,1)=tempBlock3(11:18,1)+10;
    BlockP{3} = tempBlock3;

end

% create new varianle to get probe number within block (1-10)
Block1sws(:,end+1)=Block1sws(:,2);
Block2sws(:,end+1)=Block2sws(:,2)-10;
Block3sws(:,end+1)=Block3sws(:,2)-20;
Block4sws(:,end+1)=Block4sws(:,2)-30;

[Block1SWF, Block1SWC, Block1SWB] = SW_TimeRecodingforMR (Block1sws, BlockP{1},saveWaves);
[Block2SWF, Block2SWC, Block2SWB] = SW_TimeRecodingforMR (Block2sws, BlockP{2},saveWaves);
[Block3SWF, Block3SWC, Block3SWB] = SW_TimeRecodingforMR (Block3sws, BlockP{3},saveWaves);
[Block4SWF, Block4SWC, Block4SWB] = SW_TimeRecodingforMR (Block4sws, BlockP{4},saveWaves);

% for nchan=1:max(slow_Waves(:,3))
%     [Block1SW] = SW_TimeRecodingforMR_singlechan(Block1sws, BlockP{1},saveWaves, nchan);
%     [Block2SW] = SW_TimeRecodingforMR_singlechan(Block2sws, BlockP{2},saveWaves, nchan);
%     [Block3SW] = SW_TimeRecodingforMR_singlechan(Block3sws, BlockP{3},saveWaves, nchan);
%     [Block4SW] = SW_TimeRecodingforMR_singlechan(Block4sws, BlockP{4},saveWaves, nchan);
% 
%     savefilenameB1 = ['swTimes_', subFolder '_DSWBlock1_', saveWaves.ChanLabels{nchan}, '.csv'];
%     writetable(Block1SW,[[TimingPath subFolder '/Block1'] filesep savefilenameB1]);
%     savefilenameB2 = ['swTimes_', subFolder '_DSWBlock2_', saveWaves.ChanLabels{nchan}, '.csv'];
%     writetable(Block2SW,[[TimingPath subFolder '/Block2'] filesep savefilenameB2]);
%     savefilenameB3 = ['swTimes_', subFolder '_DSWBlock3_', saveWaves.ChanLabels{nchan}, '.csv'];
%     writetable(Block3SW,[[TimingPath subFolder '/Block3'] filesep savefilenameB3]);
%     savefilenameB4 = ['swTimes_', subFolder '_DSWBlock4_', saveWaves.ChanLabels{nchan}, '.csv'];
%     writetable(Block4SW,[[TimingPath subFolder '/Block4'] filesep savefilenameB4]);
% end

%%
% SAVE

% BLOCK 1 
savefilenameB1SF = ['swTimes_', subFolder '_Block1_Fz.csv'];
writetable(Block1SWF,[[TimingPath subFolder '/Block1'] filesep savefilenameB1SF]);
savefilenameB1SC = ['swTimes_', subFolder '_Block1_Cz.csv'];
writetable(Block1SWC,[[TimingPath subFolder '/Block1'] filesep savefilenameB1SC]);
savefilenameB1SB = ['swTimes_', subFolder '_Block1_Pz.csv'];
writetable(Block1SWB,[[TimingPath subFolder '/Block1'] filesep savefilenameB1SB]);

% BLOCK 2
savefilenameB2SF = ['swTimes_', subFolder '_Block2_Fz.csv'];
writetable(Block2SWF,[[TimingPath subFolder '/Block2'] filesep savefilenameB2SF]);
savefilenameB2SC = ['swTimes_', subFolder '_Block2_Cz.csv'];
writetable(Block2SWC,[[TimingPath subFolder '/Block2'] filesep savefilenameB2SC]);
savefilenameB2SB = ['swTimes_', subFolder '_Block2_Pz.csv'];
writetable(Block2SWB,[[TimingPath subFolder '/Block2'] filesep savefilenameB2SB]);


% BLOCK 3
savefilenameB3SF = ['swTimes_', subFolder '_Block3_Fz.csv'];
writetable(Block3SWF,[[TimingPath subFolder '/Block3'] filesep savefilenameB3SF]);
savefilenameB3SC = ['swTimes_', subFolder '_Block3_Cz.csv'];
writetable(Block3SWC,[[TimingPath subFolder '/Block3'] filesep savefilenameB3SC]);
savefilenameB3SB = ['swTimes_', subFolder '_Block3_Pz.csv'];
writetable(Block3SWB,[[TimingPath subFolder '/Block3'] filesep savefilenameB3SB]);

    % BLOCK 4
savefilenameB4SF = ['swTimes_', subFolder '_Block4_Fz.csv'];
writetable(Block4SWF,[[TimingPath subFolder '/Block4'] filesep savefilenameB4SF]);
savefilenameB4SC = ['swTimes_', subFolder '_Block4_Cz.csv'];
writetable(Block4SWC,[[TimingPath subFolder '/Block4'] filesep savefilenameB4SC]);
savefilenameB4SB = ['swTimes_', subFolder '_Block4_Pz.csv'];
writetable(Block4SWB,[[TimingPath subFolder '/Block4'] filesep savefilenameB4SB]);



%end



% %%
% % Set up Event file for source reconstruction
% 
% 
% files=dir([save_path filesep 'prct_DSS_SW_*']);
% 
% for nF=1:length(files)
% 
% SubID=files(nF).name;
% SubID=regexp(SubID, '\d+', 'match');
% 
% saveWaves =load([files(nF).folder filesep files(nF).name]);
% slow_Waves= saveWaves.slow_Waves;
% 
% % Change onset and end of wave time into seconds from the probe
% Trialtimes = -25000:2:-2;
% slow_Waves(:,5)=Trialtimes(slow_Waves(:,5))/1000;
% slow_Waves(:,8)=Trialtimes(slow_Waves(:,8))/1000;
% 
% % load behaviour
% file_behav=dir([data_path filesep '..' filesep '..' filesep 'Behav' filesep 'wanderIM_behavres_s' SubID{:} '*.mat']);
% load([file_behav.folder filesep file_behav.name])
% 
% probe_res(:,end+1)=probe_res(:,5)*10-(10-probe_res(:,1)); 
% for nSW=1:length(slow_Waves)
%     thisprobenum = slow_Waves(nSW,2);
%     slow_Waves(nSW,17:21)=probe_res(probe_res(:,end)==thisprobenum,[1,5,17:19,]);
% end 
% 
% tempslow_Waves = slow_Waves(:,[2,3,5,16:21]);
% 
% FrontalElectrodes = find(ismember(saveWaves.ChanLabels,{'F1','Fz', 'F2', 'AF4', 'AF3', 'FPz', 'FP1', 'Fp2'})); %based on group level stats
% % CentralElectrodes = find(ismember(saveWaves.ChanLabels,{'FC1','FC2','CP1', 'CPz', 'CP2', 'P1','Pz','P2',}));
% % BackElectrodes = find(ismember(saveWaves.ChanLabels,{'PO3', 'POz', 'PO4','PO8','PO7', 'Oz', 'O1', 'O2'}));
% 
% %FrontalElectrodes = find(ismember(saveWaves.ChanLabels,{'Fz'})); %test 1 chan
% CentralElectrodes = find(ismember(saveWaves.ChanLabels,{'FC1','FC2','CP1', 'CPz', 'CP2', 'P1','Pz','P2',}));
% BackElectrodes = find(ismember(saveWaves.ChanLabels,{'PO3', 'POz', 'PO4','PO8','PO7', 'Oz', 'O1', 'O2'}));
% 
% Front_SWs = tempslow_Waves(ismember(tempslow_Waves(:,2),FrontalElectrodes),:);
% Center_SWs = tempslow_Waves(ismember(tempslow_Waves(:,2),CentralElectrodes),:);
% Back_SWs = tempslow_Waves(ismember(tempslow_Waves(:,2),BackElectrodes),:);
% 
% All_Front_SW_event = Front_SWs(:,[3,4]);
% All_Front_SW_event = sortrows(All_Front_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% MW_Front_SW_event = Front_SWs(Front_SWs(:,7)==2,[3,4]);
% MW_Front_SW_event = sortrows(MW_Front_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% OT_Front_SW_event = Front_SWs(Front_SWs(:,7)==1,[3,4]);
% OT_Front_SW_event = sortrows(OT_Front_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% MB_Front_SW_event = Front_SWs(Front_SWs(:,7)>=3,[3,4]);
% MB_Front_SW_event = sortrows(MB_Front_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% NOT_Front_SW_event = Front_SWs(Front_SWs(:,7)>=2,[3,4]);
% NOT_Front_SW_event = sortrows(NOT_Front_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% 
% All_Center_SW_event = Center_SWs(:,[3,4]);
% All_Center_SW_event = sortrows(All_Center_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% MW_Center_SW_event = Center_SWs(Center_SWs(:,7)==2,[3,4]);
% MW_Center_SW_event = sortrows(MW_Center_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% OT_Center_SW_event = Center_SWs(Center_SWs(:,7)==1,[3,4]);
% OT_Center_SW_event = sortrows(OT_Center_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% MB_Center_SW_event = Center_SWs(Center_SWs(:,7)>=3,[3,4]);
% MB_Center_SW_event = sortrows(MB_Center_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% NOT_Center_SW_event = Center_SWs(Center_SWs(:,7)>=2,[3,4]);
% NOT_Center_SW_event = sortrows(NOT_Center_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% 
% All_Back_SW_event = Back_SWs(:,[3,4]);
% All_Back_SW_event = sortrows(All_Back_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% MW_Back_SW_event = Back_SWs(Back_SWs(:,7)==2,[3,4]);
% MW_Back_SW_event = sortrows(MW_Back_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% OT_Back_SW_event = Back_SWs(Back_SWs(:,7)==1,[3,4]);
% OT_Back_SW_event = sortrows(OT_Back_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% MB_Back_SW_event = Back_SWs(Back_SWs(:,7)>=3,[3,4]);
% MB_Back_SW_event = sortrows(MB_Back_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% NOT_Back_SW_event = Back_SWs(Back_SWs(:,7)>=2,[3,4]);
% NOT_Back_SW_event = sortrows(NOT_Back_SW_event,[2,1],'ascend'); % here we sort rows into ascending order based on onset
% 
% 
% % Save files
% subFolder = ['s' string(SubID)];
% fprintf('Saving SW events subject: %s \n',string(SubID));
% 
% % Full directory path construction
% fullDirPath = [SWTimingPath, subFolder, filesep];
% fullDirPath = strjoin(fullDirPath, '');
% 
% % Check and create the directory if needed
% if ~exist(fullDirPath, 'dir')
%     mkdir(fullDirPath);
% end
% 
% % Variables and their names
% variables = {All_Front_SW_event, MW_Front_SW_event, OT_Front_SW_event, MB_Front_SW_event, NOT_Front_SW_event, ...
%              All_Center_SW_event, MW_Center_SW_event, OT_Center_SW_event, MB_Center_SW_event, NOT_Center_SW_event, ...
%              All_Back_SW_event, MW_Back_SW_event, OT_Back_SW_event, MB_Back_SW_event, NOT_Back_SW_event};
% variableNames = {'All_Front_SW_event', 'MW_Front_SW_event', 'OT_Front_SW_event', 'MB_Front_SW_event', 'NOT_Front_SW_event', ...
%                  'All_Center_SW_event', 'MW_Center_SW_event', 'OT_Center_SW_event', 'MB_Center_SW_event', 'NOT_Center_SW_event', ...
%                  'All_Back_SW_event', 'MW_Back_SW_event', 'OT_Back_SW_event', 'MB_Back_SW_event', 'NOT_Back_SW_event'};
% 
% variableNames = cellfun(@(x) ['events_' erase(x, '_event')], variableNames, 'UniformOutput', false);
% 
% 
% % Loop through each variable and save it
% for idx = 1:length(variables)
%     saveDataAsStruct(variables{idx}, variableNames{idx}, fullDirPath);
% end
% 
% end
% 
% 
% % Function to save variable data to a .mat file as a structure
% function saveDataAsStruct(data, fileName, fullPath)
%     % Construct the structure
%     events = struct(...
%         'label', repmat({'SW'}, size(data, 1), 1),...
%         'color', repmat({[0, 0, 0]}, size(data, 1), 1),...
%         'epochs', num2cell(data(:,2)),...
%         'times', num2cell(data(:,1))...
%     );
%     
%     % Save the structure to a .mat file
%     save(fullfile(fullPath, [fileName, '.mat']), 'events');
% end
% 
% 

% 
% % Combine all blocks into cell arrays
% blocks_F = {Block1SWF, Block2SWF, Block3SWF, Block4SWF};
% blocks_C = {Block1SWC, Block2SWC, Block3SWC, Block4SWC};
% blocks_B = {Block1SWB, Block2SWB, Block3SWB, Block4SWB};
% 
% % Initialize vectors for the plot
% x_F = [];
% y_F = [];
% x_C = [];
% y_C = [];
% x_B = [];
% y_B = [];
% 
% currentOffset = 0;
% breakDuration = 00;
% % Construct plot data for blocks_F
% for b = 1:length(blocks_F)
%     block = blocks_F{b};
%     onsets = block.Onset;
%     durations = block.Duration;
%     amplitudes = block.Amplitude;
%     
%     for i = 1:height(block)
%         x_F = [x_F, currentOffset + onsets(i), currentOffset + onsets(i), currentOffset + onsets(i) + durations(i), currentOffset + onsets(i) + durations(i)];
%         y_F = [y_F, 0, amplitudes(i), amplitudes(i), 0];
%     end
%     
%     currentOffset = currentOffset + max(onsets + durations) + breakDuration;
% 
% end
% x_F = [0, x_F, max(x_F)];
% y_F = [0, y_F, 0];
% 
% % Reset offset for the next blocks
% currentOffset = 0;
% 
% % Construct plot data for blocks_C
% for b = 1:length(blocks_C)
%     block = blocks_C{b};
%     onsets = block.Onset;
%     durations = block.Duration;
%     amplitudes = block.Amplitude;
%     
%     for i = 1:height(block)
%         x_C = [x_C, currentOffset + onsets(i), currentOffset + onsets(i), currentOffset + onsets(i) + durations(i), currentOffset + onsets(i) + durations(i)];
%         y_C = [y_C, 0, amplitudes(i), amplitudes(i), 0];
%     end
%     
%      currentOffset = currentOffset + max(onsets + durations) + breakDuration;
% end
% 
% x_C = [0, x_C, max(x_C)];
% y_C = [0, y_C, 0];
% 
% % Reset offset for the next blocks
% currentOffset = 0;
% 
% % Construct plot data for blocks_B
% for b = 1:length(blocks_B)
%     block = blocks_B{b};
%     onsets = block.Onset;
%     durations = block.Duration;
%     amplitudes = block.Amplitude;
%     
%     for i = 1:height(block)
%         x_B = [x_B, currentOffset + onsets(i), currentOffset + onsets(i), currentOffset + onsets(i) + durations(i), currentOffset + onsets(i) + durations(i)];
%         y_B = [y_B, 0, amplitudes(i), amplitudes(i), 0];
%     end
%     
%     currentOffset = currentOffset + max(onsets + durations) + breakDuration;
% end
% x_B = [0, x_B, max(x_B)];
% y_B = [0, y_B, 0];
% 
% 
% % Plot the data
% figure;
% hold on;
% plot(x_F, y_F, 'LineWidth', 2, 'Color', [0.7, 0.4, 0.4, 0.5]); % Blue for Block*SWF
% %plot(x_C, y_C, 'LineWidth', 2, 'Color', [0.4, 0.7, 0.4, 0.5]); % Red for Block*SWC
% %plot(x_B, y_B, 'LineWidth', 2, 'Color', [0.4, 0.4, 0.7, 0.5]); % Green for Block*SWB
% xlabel('Time (s)');
% ylabel('Amplitude');
% title(['Regressor Plot Across Blocks - ',subFolder]);
% %legend('Block*SWF', 'Block*SWC', 'Block*SWB');
% grid on;
% hold off;

end