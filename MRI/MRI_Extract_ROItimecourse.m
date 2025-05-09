%%

% marsbar raw data
% mpcorrected are eigen values
% whitefilt has big jumps across blocks



%%

% TimeCourseDataset_ns = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/b*_timecourse_n*.mat'); % fr*_timecourse_s*.mat 
TimeCourseDataset = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/b*_timecourse_s*.mat'); % fr*_timecourse_s*.mat 

% TimeCourseDataset = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/ot*.mat');
% TimeCourseDataset = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Old/rffa_timecourse_adjunsted*.mat');
% TimeCourseDataset = TimeCourseDataset([1:25],:);
% 
% TimeCourseDataset_B1 = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/VOI_f*mp*_1.mat'); % VOI_f*whitefilt*_1 mpcorrected
% TimeCourseDataset_B2 = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/VOI_f*mp*_2.mat');
% TimeCourseDataset_B3 = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/VOI_f*mp*_3.mat');
% TimeCourseDataset_B4 = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/VOI_f*mp*_4.mat');

%GrandaverageEpoch = [];
for k=1:size(TimeCourseDataset,1)

    TCpathparts = split(TimeCourseDataset(k).name, "_");
    %TCpathparts = split(TimeCourseDataset_B1(k).name, "_");
    % 
    % if length(TCpathparts)>5 && length(TCpathparts)<7
    %     filen = [TCpathparts{1}, '_', TCpathparts{2}, '_', TCpathparts{3}];
    %      SubID = TCpathparts{5};
    % elseif length(TCpathparts)>6 && length(TCpathparts)<8
    %     %filen = [TCpathparts{2}, '_', TCpathparts{3}, '_', TCpathparts{7}];
    %     filen = [TCpathparts{1}, '_', TCpathparts{2}];
    %     SubID = TCpathparts{7};
    % elseif length(TCpathparts)>=8
    %     %filen = [TCpathparts{2}, '_', TCpathparts{3}, '_', TCpathparts{7}];
    %     filen = [TCpathparts{2}, '_', TCpathparts{3} ];
    %     SubID = TCpathparts{6};
    % else 
    %     filen = [TCpathparts{1}, '_', TCpathparts{2}];
    %     SubID = TCpathparts{5};
    % end

    filen = [TCpathparts{1}, '_', TCpathparts{2} ];
    SubID = TCpathparts{5};

   
    load([TimeCourseDataset(k).folder filesep TimeCourseDataset(k).name]);
    TimeCourseM = ROI_data_mars;
    % 
    % load([TimeCourseDataset_ns(k).folder filesep TimeCourseDataset_ns(k).name]);
    % TimeCourseM_ns = ROI_data_mars;
    
    % TimeCourse = [];
    % load([TimeCourseDataset_B1(k).folder filesep TimeCourseDataset_B1(k).name]);
    % TimeCourse = [TimeCourse; Y];
    % load([TimeCourseDataset_B2(k).folder filesep TimeCourseDataset_B2(k).name]);
    % TimeCourse = [TimeCourse; Y];
    % load([TimeCourseDataset_B3(k).folder filesep TimeCourseDataset_B3(k).name]);
    % TimeCourse = [TimeCourse; Y];
    % load([TimeCourseDataset_B4(k).folder filesep TimeCourseDataset_B4(k).name]);
    % TimeCourse = [TimeCourse; Y];
    % TimeCourse = TimeCourse';

    TimeCourse = TimeCourseM;


    BlockDataset = dir(['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/TimingforfMRI/',SubID, '/*/full*.csv']); 
    
    Block1 = readtable([BlockDataset(1).folder filesep BlockDataset(1).name]);
    Block2 = readtable([BlockDataset(2).folder filesep BlockDataset(2).name]);
    Block3 = readtable([BlockDataset(3).folder filesep BlockDataset(3).name]);
    Block4 = readtable([BlockDataset(4).folder filesep BlockDataset(4).name]);
    
    Block2.Onset = Block2.Onset+Block1.Onset(end,1);
    Block3.Onset = Block3.Onset+Block2.Onset(end,1);
    Block4.Onset = Block4.Onset+Block3.Onset(end,1);
    AllBlocks = [Block1; Block2; Block3; Block4];


%     % Identify rows where Type starts with 'ProbeEnd'
%     probeEndRows = startsWith(AllBlocks.Type, 'ProbeEnd');
%     % Find the indices of these rows
%     probeEndIndices = find(probeEndRows);
%     % Initialize an array to store indices of rows to be included in the new table
%     followIndices = [];
%     % Loop through each identified ProbeEnd row to get the following row
%     for i = 1:length(probeEndIndices)
%         % Check if the current ProbeEnd row is not the last row
%         if probeEndIndices(i) < height(AllBlocks)
%             % Add the index of the row that follows the current ProbeEnd row
%             tempInd = probeEndIndices(i) + 1;
%             followIndices = [followIndices; tempInd'];
%         end
%     end
%     % Create a new table with the selected rows
%     FirstFaceTrials = AllBlocks(followIndices, :);


%     SWDataset = dir(['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/TimingforfMRI/',SubID, '/*/*_Block*Cz.csv']); 
%     SWpathparts = split(SWDataset(1).name, "_");
%     ChanID = split(SWpathparts{4}, ".");
%     ChanID = ChanID{1};
% 
% 
%     SWBlock1 = readtable([SWDataset(1).folder filesep SWDataset(1).name]);
%     SWBlock2 = readtable([SWDataset(2).folder filesep SWDataset(2).name]);
%     SWBlock3 = readtable([SWDataset(3).folder filesep SWDataset(3).name]);
%     SWBlock4 = readtable([SWDataset(4).folder filesep SWDataset(4).name]);
%     
%     SWBlock2.Onset = SWBlock2.Onset+SWBlock1.Onset(end,1);
%     SWBlock3.Onset = SWBlock3.Onset+SWBlock2.Onset(end,1);
%     SWBlock4.Onset = SWBlock4.Onset+SWBlock3.Onset(end,1);
%     AllSWBlocks = [SWBlock1; SWBlock2; SWBlock3; SWBlock4];



%% FILTER
DataToFilt = TimeCourse;
% Sample rate 
Fs = 1 / 1.56;  % Sampling frequency (Hz)
% Filter design
lowcut = 0.01;  % Lower frequency bound (Hz)
highcut = 0.1;  % Upper frequency bound (Hz)
nyquist = Fs/2; % Nyquist frequency
% Calculate normalized frequency for filter design
Wn = [lowcut highcut] / nyquist;
% Design a Butterworth bandpass filter
% (Change the order as needed, here we use 4th order)
[b, a] = butter(4, Wn, 'bandpass');
% Apply the filter using filtfilt for zero-phase filtering
filteredData = filtfilt(b, a, DataToFilt);

%TimeCourse = filteredData;
TimeCourse = TimeCourse - mean(TimeCourse,"all");

    % Select probe onset
    ProbeStop = {'ProbeStopONTASK', 'ProbeStopMINDWANDER', 'ProbeStopMINDBLANK', 'ProbeStopNORECALL'};
    isProbeStop = ismember(AllBlocks.Type, ProbeStop);
    SubsetTimes = AllBlocks(isProbeStop,:);

    TR = 1.56; % Repetition Time in seconds
    epochWindow = 20; % Time window in seconds
    epochPoints = round(epochWindow / TR); % Number of points in the time window

    % Extract event onsets
    Tagettable = SubsetTimes; %AllSWBlocks;
    onsets = Tagettable.Onset;   
    numEvents = height(Tagettable);
    % Initialize matrix to store all epochs
    allEpochs = NaN(numEvents, epochPoints + 1); % NaN(numEvents, 2 * epochPoints + 1);
    
    for i = 1:numEvents
        % Convert onset time to data point index
        eventIndex = round(onsets(i) / TR);

        
        % Define start and end points of the epoch
        startIndex = eventIndex - epochPoints;
        endIndex = eventIndex; %eventIndex + epochPoints;
        
        % Extract the epoch, using NaNs for out-of-bounds indices
        epoch = NaN(1, epochPoints + 1); %NaN(1, 2 * epochPoints + 1);
        validStart = max(1, startIndex);
        validEnd = min(length(TimeCourse), endIndex);
        epoch((validStart - startIndex + 1):(validEnd - startIndex + 1)) = TimeCourse(validStart:validEnd);
        
        % Store the epoch in the matrix
        allEpochs(i, :) = epoch;
        
    end

    % save 
    MeanEpochs = mean(allEpochs,2);
    ROIsavepath ='/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI';
    savefilename1 = [ROIsavepath filesep 'Preprobe20sARV_',  filen, '_mars2_', SubID '.mat'];
    save(savefilename1,'MeanEpochs')
    fprintf('Done subject: %s \n',SubID);


%     % remove all mean activity in epoch
%     allEpochs(:,:)= allEpochs(:,:)-mean(TimeCourse);
%     %allEpochs(:,:)= allEpochs(:,:)-allEpochs(:,14);
% 
%     % Calculate the average activity across all epochs
%     averageEpoch = nanmean(allEpochs, 1);
% 
%     GrandaverageEpoch(k, :) = averageEpoch;
%      % Calculate the standard error of the mean
%     sem = nanstd(allEpochs, 0, 1) / sqrt(numEvents);
%     % Calculate the 95% confidence interval
%     ci95 = sem * 1.96;
% 
%     % Plot all individual epochs
%     figure;
%     hold on;
%     timeAxis = linspace(-epochWindow, epochWindow, size(allEpochs, 2));
%     for i = 1:numEvents
%         plot(timeAxis, allEpochs(i, :), 'Color', [0.8, 0.8, 0.8]); % Light gray color for individual epochs
%     end
%     
%     % Plot the average activity
%     plot(timeAxis, averageEpoch, 'b', 'LineWidth', 2); % Blue line for average
%     
%     % Plot the 95% confidence interval
%     fill([timeAxis, fliplr(timeAxis)], [averageEpoch + ci95, fliplr(averageEpoch - ci95)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     
%     title(['SW-locked Average Activity' newline SubID ' - ' TCpathparts{1}, ' ', TCpathparts{2}, ' based on  ', ChanID]);
%     xlabel('Time course [0 = SW onset] (sec)')
%     ylabel('BOLD activity')
%     set(gcf,'Color','w')
%     set(gca,'FontSize',18,'FontWeight','bold')
%     grid on;
%     hold off;
%     saveas(gcf,['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/ROIFigures/' filesep SubID '_' TCpathparts{1}, '_', TCpathparts{2}, '_', ChanID, '.png'])
% 

%     figure
%     xax = (0:length(TimeCourse)-1) * 1.56;
%     plot(xax, TimeCourse);
%     hold on
%     title(['Full time course' newline SubID ' - ' TCpathparts{1}, ' ', TCpathparts{2}, ' based on  ', ChanID]);
%     xlabel('Time course (sec)')
%     ylabel('Average BOLD activity')
%     axis(gca, 'tight')  
% 
%         % Add vertical lines for trial onsets
%     for i = 1:height(AllSWBlocks)
%         x = AllSWBlocks.Onset(i);
%         line([x, x], ylim, 'Color', [0.4, 0.4, 0.4, 0.5], 'LineStyle', '-', 'LineWidth', 0.5);
%     end
% 
%     saveas(gcf,['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/ROIFigures' filesep SubID '_TimeCoursewithSW_', ChanID, '.png'])


end

% 
%  Grandaverage = mean(GrandaverageEpoch, 1);
%      % Calculate the standard error of the mean
%  GAsem = std(GrandaverageEpoch, 0, 1) / sqrt(length(TimeCourseDataset));
%     % Calculate the 95% confidence interval
%  GAci95 = GAsem * 1.96;
% 
%  % Plot all individual epochs
% figure;
% hold on;
% timeAxis = linspace(-epochWindow, epochWindow, size(Grandaverage, 2));
% 
% % Plot the average activity
% plot(timeAxis, Grandaverage, 'b', 'LineWidth', 2); % Blue line for average
% 
% % Plot the 95% confidence interval
% fill([timeAxis, fliplr(timeAxis)], [Grandaverage + GAci95, fliplr(Grandaverage - GAci95)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
% title(['SW-locked Grand Average Activity' newline TCpathparts{1}, ' ', TCpathparts{2}, ' based on  ', ChanID]);
% xlabel('Time course [0 = SW onset] (sec)')
% ylabel('BOLD activity')
% set(gcf,'Color','w')
% set(gca,'FontSize',18,'FontWeight','bold')
% grid on;
% hold off;
% saveas(gcf,['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Figures/ROIFigures' filesep 'GrandAverage_' TCpathparts{1}, '_', TCpathparts{2}, '_', ChanID, '.png'])
% 

%%
% figure;
% plot(TimeCourseM)
% hold on
% plot(TimeCourse)
% plot(filteredData)

% TimeCourseMnorm = TimeCourseM - mean(TimeCourseM,"all");
% TimeCourseMnorm_ns = TimeCourseM_ns - mean(TimeCourseM_ns,"all");
% 
% 
% figure;
% plot(TimeCourseMnorm, 'Color', [0.5 0.0 0.5 0.4], 'LineWidth', 2)
% hold on
% plot(TimeCourseMnorm_ns, 'Color', [0 0.5 0.5 0.4],'LineWidth', 2)
% plot(filteredData, 'Color', [0 0.5 0 0.4],'LineWidth', 2)
% plot(filteredData_ns, 'Color', [0.7 0.5 0 0.4],'LineWidth', 2)
% 
% 
% xlabel('Time course'); % Label for the x-axis
% ylabel('Bold activity'); % Label for the y-axis
% %legend('Marsbar', 'SPM');
% legend('Raw timecourse', 'Raw timecourse - non-smoothing', 'Filtered timecourse', 'Filtered timecourse - non-smoothing');
% hold off; % Release the hold
% format_fig;