



%%

TimeCourseDataset = dir('/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/ROI/Raw time course/f*_L*.mat');

GrandaverageEpoch = [];
for k=1:size(TimeCourseDataset,1)

    TCpathparts = split(TimeCourseDataset(k).name, "_");
    SubID = split(TCpathparts{4}, ".");
    SubID = SubID{1};

    load([TimeCourseDataset(k).folder filesep TimeCourseDataset(k).name]);

    TimeCourse = ROI_data_mars;

    BlockDataset = dir(['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/TimingforfMRI/',SubID, '/*/full*.csv']); 
    
    Block1 = readtable([BlockDataset(1).folder filesep BlockDataset(1).name]);
    Block2 = readtable([BlockDataset(2).folder filesep BlockDataset(2).name]);
    Block3 = readtable([BlockDataset(3).folder filesep BlockDataset(3).name]);
    Block4 = readtable([BlockDataset(4).folder filesep BlockDataset(4).name]);
    
    Block2.Onset = Block2.Onset+Block1.Onset(end,1);
    Block3.Onset = Block3.Onset+Block2.Onset(end,1);
    Block4.Onset = Block4.Onset+Block3.Onset(end,1);
    AllBlocks = [Block1; Block2; Block3; Block4];

    SWDataset = dir(['/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/TimingforfMRI/',SubID, '/*/*_Block*Pz.csv']); 
    SWpathparts = split(SWDataset(1).name, "_");
    ChanID = split(SWpathparts{4}, ".");
    ChanID = ChanID{1};


    SWBlock1 = readtable([SWDataset(1).folder filesep SWDataset(1).name]);
    SWBlock2 = readtable([SWDataset(2).folder filesep SWDataset(2).name]);
    SWBlock3 = readtable([SWDataset(3).folder filesep SWDataset(3).name]);
    SWBlock4 = readtable([SWDataset(4).folder filesep SWDataset(4).name]);
    
    SWBlock2.Onset = SWBlock2.Onset+SWBlock1.Onset(end,1);
    SWBlock3.Onset = SWBlock3.Onset+SWBlock2.Onset(end,1);
    SWBlock4.Onset = SWBlock4.Onset+SWBlock3.Onset(end,1);
    AllSWBlocks = [SWBlock1; SWBlock2; SWBlock3; SWBlock4];

    TRtimes = 1:size(TimeCourse,1);
    TRtimes = TRtimes*1.56;

    TR = 1.56; % Repetition Time in seconds
    epochWindow = 20; % Time window in seconds
    epochPoints = round(epochWindow / TR); % Number of points in the time window

    % Extract event onsets
    onsets = AllSWBlocks.Onset;
    numEvents = height(AllSWBlocks);
    allEpochs = [];
    
    for i = 1:numEvents
        % Convert onset time to data point index
        eventIndex = round(onsets(i) / TR);
        
        % Define start and end points of the epoch
        startIndex = eventIndex - epochPoints;
        endIndex = eventIndex + epochPoints;
        
        % Ensure indices are within bounds
        if startIndex < 1
            startIndex = 1;
        endIndex = startIndex + 2 * epochPoints; % Maintain epoch length
        end
        if endIndex > length(TimeCourse)
            endIndex = length(TimeCourse);
            startIndex = endIndex - 2 * epochPoints; % Maintain epoch length
        end
        
        % Extract the epoch
        allEpochs(i, :) = TimeCourse(startIndex:endIndex);
    end

    %removing probe onset activity
    allEpochs(:,:)= allEpochs(:,:)-allEpochs(14,:);
    % remove all mean activity
    %allEpochs(:,:)= allEpochs(:,:)-mean(allEpochs(:,:));

    % Calculate the average activity across all epochs
    averageEpoch = mean(allEpochs, 1);

    GrandaverageEpoch(k, :) = averageEpoch;
     % Calculate the standard error of the mean
    sem = std(allEpochs, 0, 1) / sqrt(numEvents);
    % Calculate the 95% confidence interval
    ci95 = sem * 1.96;



    % Plot all individual epochs
    figure;
    hold on;
    timeAxis = linspace(-epochWindow, epochWindow, size(allEpochs, 2));
    for i = 1:numEvents
        plot(timeAxis, allEpochs(i, :), 'Color', [0.8, 0.8, 0.8]); % Light gray color for individual epochs
    end
    
    % Plot the average activity
    plot(timeAxis, averageEpoch, 'b', 'LineWidth', 2); % Blue line for average
    
    % Plot the 95% confidence interval
    fill([timeAxis, fliplr(timeAxis)], [averageEpoch + ci95, fliplr(averageEpoch - ci95)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    title(['SW-locked Average Activity' newline SubID ' - ' TCpathparts{1}, ' ', TCpathparts{2}, ' based on  ', ChanID]);
    xlabel('Time course [0 = SW onset] (sec)')
    ylabel('BOLD activity')
    set(gcf,'Color','w')
    set(gca,'FontSize',18,'FontWeight','bold')
    grid on;
    hold off;
    saveas(gcf,['/Users/kuszti/Documents/MATLAB/EEG_FMRI/Figures/ROIFigures' filesep SubID '_' TCpathparts{1}, '_', TCpathparts{2}, '_', ChanID, '.png'])

end


 Grandaverage = mean(GrandaverageEpoch, 1);
     % Calculate the standard error of the mean
 GAsem = std(GrandaverageEpoch, 0, 1) / sqrt(length(TimeCourseDataset));
    % Calculate the 95% confidence interval
 GAci95 = GAsem * 1.96;

 % Plot all individual epochs
figure;
hold on;
timeAxis = linspace(-epochWindow, epochWindow, size(Grandaverage, 2));

% Plot the average activity
plot(timeAxis, Grandaverage, 'b', 'LineWidth', 2); % Blue line for average

% Plot the 95% confidence interval
fill([timeAxis, fliplr(timeAxis)], [Grandaverage + GAci95, fliplr(Grandaverage - GAci95)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

title(['SW-locked Grand Average Activity' newline TCpathparts{1}, ' ', TCpathparts{2}, ' based on  ', ChanID]);
xlabel('Time course [0 = SW onset] (sec)')
ylabel('BOLD activity')
set(gcf,'Color','w')
set(gca,'FontSize',18,'FontWeight','bold')
grid on;
hold off;

saveas(gcf,['/Users/kuszti/Documents/MATLAB/EEG_FMRI/Figures/ROIFigures' filesep 'GrandAverage_' TCpathparts{1}, '_', TCpathparts{2}, '_', ChanID, '.png'])
