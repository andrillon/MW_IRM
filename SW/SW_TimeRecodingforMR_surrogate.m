function [SimFTable, SimCTable, SimBTable] = SW_TimeRecodingforMR_surrogate (BlockSW, BlockP, saveWaves) 

% set parameters to extract
onsetp = 5; % 5=start neg. wave; 8=peak neg wave;6=start pos wave; 10=peak pos wave
durp = 7; %7=end of wave
ampp = 4; % 9=neg peak, 11=pos peak, 4=peak2peak


% % create a tabel on onsets,durations,amplitudes
OnsetT = table;
DurT = table;
AmpT = table;
SWProbeT = table;

nPlist=unique(BlockSW(:,end))';
for nP=nPlist % for every probe 
    
    allOnset ={};
    allDuration ={};
    allAmplitude ={};

    for nEl=unique(BlockSW(:,3))' %and every electrode 
        onset=[];
        duration =[];
        amplitude = [];
        temp_waves=BlockSW(BlockSW(:,end)==nP & BlockSW(:,3)==nEl,:);% take all slow waves
        for nW=1:size(temp_waves,1)
            onset = [onset; BlockP(BlockP(:,1)==nP,2) + (temp_waves(nW,onsetp)/1000)];%subtract negative onset value from probe onset in seconds
            duration = [duration; (temp_waves(nW,durp)-temp_waves(nW,onsetp))/1000]; % calculate duration in seconds
            amplitude = [amplitude; abs(temp_waves(nW,ampp))]; % take absolute amplitude of negative peak
            
        end
        allOnset{nEl}=onset; % for each electrode, we then copy all onset times 
        allDuration{nEl} =duration;
        allAmplitude{nEl} =amplitude;
    end

    max_len = max(cellfun(@length,allOnset)); % get the longest vector length

    for k = 1:size(allOnset,2)
        allOnset{k} = [allOnset{k}; NaN(max_len-size(allOnset{k},1),1)]; % fill each vector with missing to have the longest vector length
        allDuration{k} = [allDuration{k}; NaN(max_len-size(allDuration{k},1),1)]; % fill each vector with missing to have the longest vector length
        allAmplitude{k} = [allAmplitude{k}; NaN(max_len-size(allAmplitude{k},1),1)]; % fill each vector with missing to have the longest vector length
    end
OnsetT = [OnsetT; array2table(cell2mat(allOnset))]; 
DurT = [DurT; array2table(cell2mat(allDuration))]; 
AmpT = [AmpT; array2table(cell2mat(allAmplitude))]; 

% nnz(~isnan(table2array(OnsetT))) to check if number of SWs match BLOCKSW

% record corresponding probe responses for each onset
SWProbeT = [SWProbeT; array2table(repmat(BlockP(BlockP(:,1)==nP,[1,5:7]),size(array2table(cell2mat(allOnset)),1),1))]; 
end

% outcome matrix: 63 Electrodes all slow waves ( inlcudes NAN rows)
OnsetT.Properties.VariableNames=saveWaves.ChanLabels;
DurT.Properties.VariableNames=saveWaves.ChanLabels;
AmpT.Properties.VariableNames=saveWaves.ChanLabels;
SWProbeT.Properties.VariableNames={'ProbeNum', 'Q1', 'Q2', 'Q3'};


% % select subset of electrodes
% FrontalElectrodes = ismember(saveWaves.ChanLabels,{'FPz', 'FP1', 'Fp2'});
% CentralElectrodes = ismember(saveWaves.ChanLabels,{'C1','Cz', 'C2'});
% BackElectrodes = ismember(saveWaves.ChanLabels,{'Oz', 'O1', 'O2'});

FrontalElectrodes = ismember(saveWaves.ChanLabels,{'FPz'});
CentralElectrodes = ismember(saveWaves.ChanLabels,{'Cz'});
BackElectrodes = ismember(saveWaves.ChanLabels,{'Oz'});

% FrontalElectrodes = ismember(saveWaves.ChanLabels,{'F1','Fz', 'F2', 'AF4', 'AF3', 'FPz', 'FP1', 'Fp2'}); %based on group level stats
% CentralElectrodes = ismember(saveWaves.ChanLabels,{'FC1','FC2','CP1', 'CPz', 'CP2', 'P1','Pz','P2',});
% BackElectrodes = ismember(saveWaves.ChanLabels,{'PO3', 'POz', 'PO4','PO8','PO7', 'Oz', 'O1', 'O2'});
% 
% FrontalElectrodes = ismember(saveWaves.ChanLabels,{'AF7','Fz', 'AF8', 'AF4', 'AF3', 'FPz', 'FP1', 'Fp2'}); %based on group level stats
% CentralElectrodes = ismember(saveWaves.ChanLabels,{'FC1','FC2','CP1', 'CPz', 'CP2', 'C1', 'Cz', 'C2', 'P1','Pz','P2',});
% BackElectrodes = ismember(saveWaves.ChanLabels,{'PO3', 'POz', 'PO4','PO8','PO7', 'Oz', 'O1', 'O2'});
% 
% VERSION A: inlcude all waves

% % we need to get timecourse 
vOnset_F = reshape(table2array(OnsetT(:,FrontalElectrodes)), 1, [])';
vOnset_F(isnan(vOnset_F))=[];
vDur_F = reshape(table2array(DurT(:,FrontalElectrodes)), 1, [])';
vDur_F(isnan(vDur_F))=[];
vAmp_F= reshape(table2array(AmpT(:,FrontalElectrodes)), 1, [])';
vAmp_F(isnan(vAmp_F))=[];
% keep reshaping the probe table 
vSWP_F =[];
vSWP_F = repmat(table2array(SWProbeT), sum(FrontalElectrodes),1); 
tempOnset_F = reshape(table2array(OnsetT(:,FrontalElectrodes)), 1, [])';
vSWP_F(isnan(tempOnset_F),:)=[];

temp_F = [vOnset_F, vDur_F, vAmp_F, vSWP_F];
temp_F = sortrows(temp_F,[1,3],'ascend');

vOnset_C = reshape(table2array(OnsetT(:,CentralElectrodes)), 1, [])';
vOnset_C(isnan(vOnset_C))=[];
vDur_C = reshape(table2array(DurT(:,CentralElectrodes)), 1, [])';
vDur_C(isnan(vDur_C))=[];
vAmp_C= reshape(table2array(AmpT(:,CentralElectrodes)), 1, [])';
vAmp_C(isnan(vAmp_C))=[];
% keep reshaping the probe table 
vSWP_C =[];
vSWP_C = repmat(table2array(SWProbeT), sum(CentralElectrodes),1); 
tempOnset_C = reshape(table2array(OnsetT(:,CentralElectrodes)), 1, [])';
vSWP_C(isnan(tempOnset_C),:)=[];

temp_C = [vOnset_C, vDur_C, vAmp_C, vSWP_C];
temp_C = sortrows(temp_C,[1,3],'ascend');


vOnset_B = reshape(table2array(OnsetT(:,BackElectrodes)), 1, [])';
vOnset_B(isnan(vOnset_B))=[];
vDur_B = reshape(table2array(DurT(:,BackElectrodes)), 1, [])';
vDur_B(isnan(vDur_B))=[];
vAmp_B= reshape(table2array(AmpT(:,BackElectrodes)), 1, [])';
vAmp_B(isnan(vAmp_B))=[];
% keep reshaping the probe table 
vSWP_B =[];
vSWP_B = repmat(table2array(SWProbeT), sum(BackElectrodes),1); 
tempOnset_B = reshape(table2array(OnsetT(:,BackElectrodes)), 1, [])';
vSWP_B(isnan(tempOnset_B),:)=[];

temp_B = [vOnset_B, vDur_B, vAmp_B, vSWP_B];
temp_B = sortrows(temp_B,[1,3],'ascend'); % here we sort rows into ascending order based on onset and secondly based on amplitude

% Need better way but for now removing duplicate onsets (smaller amplitudes
% are removed
dx1 = find(diff(temp_F(:,1))==0);
temp_F(dx1,:)=[];
dx2 = find(diff(temp_C(:,1))==0);
temp_C(dx2,:)=[];
dx3 = find(diff(temp_B(:,1))==0);
temp_B(dx3,:)=[];

FTable = array2table(temp_F, 'VariableNames', {'Onset', 'Duration', 'Amplitude', 'ProbeNum', 'Q1', 'Q2', 'Q3'});
CTable = array2table(temp_C, 'VariableNames', {'Onset', 'Duration', 'Amplitude', 'ProbeNum', 'Q1', 'Q2', 'Q3'});
BTable = array2table(temp_B, 'VariableNames', {'Onset', 'Duration', 'Amplitude', 'ProbeNum', 'Q1', 'Q2', 'Q3'});


%% Simlate onset, duration and amplidute values based on selected data

SWTables{1} = FTable;
SWTables{2} = CTable;
SWTables{3} = BTable;

for j = 1:size(SWTables,2)

    targetTable = SWTables{j};
    % Initialize the simulated dataset table
    numRows = height(targetTable);
    SimulatedTable = table(NaN(numRows, 1), NaN(numRows, 1), NaN(numRows, 1), 'VariableNames', {'Onset', 'Duration', 'Amplitude'});
    window_before_probes=20; % in seconds
    % Compute the mean and standard deviation for Duration and Amplitude
    meanDuration = mean(targetTable.Duration);
    stdDuration = std(targetTable.Duration);
    meanAmplitude = mean(targetTable.Amplitude);
    stdAmplitude = std(targetTable.Amplitude);
    temp_sim = NaN(numRows, 3);
    
    
    for i = 1:numRows
        % Get the probe onset times and calculate the allowed onset window
        probeOnset = BlockP(:, 2)';
        allowedWindowStart = probeOnset - window_before_probes;
        allowedWindowEnd = probeOnset;
        
        % Remove negative start times
        allowedWindowStart(allowedWindowStart < 0) = 0;
    
        % Ensure no overlap with existing data
        overlap = true;
        while overlap
            % Generate a random onset time within the allowed window
            valid = false;
            while ~valid
                simonset = allowedWindowStart(randi(length(allowedWindowStart))) + rand * (window_before_probes);
                % Ensure onset is within the valid range
                if any(simonset >= allowedWindowStart & simonset <= allowedWindowEnd)
                    valid = true;
                end
            end
    
            valid2 = false;
            while ~valid2
                % Generate a random duration and amplitude based on the mean and std
                simduration = meanDuration + stdDuration * randn;
                simamplitude = meanAmplitude + stdAmplitude * randn;
                % check that duration and amplitude is larger than 0
                if simduration > 0 && simamplitude > 0
                valid2 = true;
                end
            end
    

        % Check for overlap with existing events in targetTable
        overlap_with_existing = any((simonset >= targetTable.Onset & simonset <= targetTable.Onset + targetTable.Duration) | ...
                                    (simonset + simduration >= targetTable.Onset & simonset + simduration <= targetTable.Onset + targetTable.Duration) | ...
                                    (simonset <= targetTable.Onset & simonset + simduration >= targetTable.Onset + targetTable.Duration));

        % Check for overlap with already generated simulated events
        overlap_with_simulated = any((simonset >= temp_sim(:, 1) & simonset <= temp_sim(:, 1) + temp_sim(:, 2)) | ...
                                     (simonset + simduration >= temp_sim(:, 1) & simonset + simduration <= temp_sim(:, 1) + temp_sim(:, 2)) | ...
                                     (simonset <= temp_sim(:, 1) & simonset + simduration >= temp_sim(:, 1) + temp_sim(:, 2)));


        % Update the overlap
        overlap = overlap_with_existing || overlap_with_simulated;

        end
    
        temp_sim(i,1) = simonset;
        temp_sim(i,2) = simduration;
        temp_sim(i,3) = simamplitude;
        
    end
    
        temp_sim = sortrows(temp_sim,1,'ascend'); % here we sort rows into ascending order based on onset
   
         % Store the simulated data
        SimulatedTable.Onset = temp_sim(:,1);
        SimulatedTable.Duration = temp_sim(:,2);
        SimulatedTable.Amplitude = temp_sim(:,3);
        SimSWTable{j} = SimulatedTable;
end

SimFTable = SimSWTable{1};
SimCTable = SimSWTable{2};
SimBTable = SimSWTable{3};

end