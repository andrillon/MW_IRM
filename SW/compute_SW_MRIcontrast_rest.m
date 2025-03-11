 
%% Init

clear all;
close all;
cd /Users/kuszti/Documents/Git/MW_IRM/SW/
run ../localdef.m

files=dir([save_path filesep 'prct_DSS_RS_SW_M*']);
EEGfiles=dir([data_path filesep filesep 'RestingState/*clean_nt.set']);

for nF=1:length(files)

saveWaves =load([files(nF).folder filesep files(nF).name]);
slow_Waves= saveWaves.slow_Waves;

SubID=files(nF).name;
SubID=SubID(21:end-4);
subFolder = ['s' (SubID)];
fprintf('Running subject: %s \n',subFolder);

 % load EEG
 EEG = pop_loadset( 'filename',[EEGfiles(nF).folder filesep EEGfiles(nF).name]);

 uEvent = EEG.urevent;
 % get epoched event file
 Event = struct2table(EEG.event);

% set parameters to extract
onsetp = 5; % 5=start neg. wave; 8=peak neg wave;6=start pos wave; 10=peak pos wave
durp = 7; %7=end of wave
ampp = 4; % 9=neg peak, 11=pos peak, 4=peak2peak

% % create a tabel on onsets,durations,amplitudes
OnsetT = table;
DurT = table;
AmpT = table;
nPlist=unique(slow_Waves(:,16))';
for nP=nPlist % for every probe 
    
    allOnset ={};
    allDuration ={};
    allAmplitude ={};

    %select epoch events
    tempEvent = Event(Event.epoch(:,:)==nP,:);
    % the EEG is epoched from the first S1 event in 30sec long consequitive
    % epochs so the first X event in uEvent is the first S1 event 
    epochstartInd = table2array(tempEvent(strcmp('X', tempEvent.type),9));
    epochstartLatency = uEvent(epochstartInd).latency;

    for nEl=unique(slow_Waves(:,3))' %and every electrode 
        onset=[];
        duration =[];
        amplitude = [];
        temp_waves=slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,3)==nEl,:);% take all slow waves
        for nW=1:size(temp_waves,1)
            onset = [onset; epochstartLatency + (temp_waves(nW,onsetp))];% add onset values to time difference from first scaner in latency
            duration = [duration; (temp_waves(nW,durp)-temp_waves(nW,onsetp))]; % calculate duration
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

% nnz(~isnan(table2array(OnsetT))) to check if number of SWs match 
end

% outcome matrix: 63 Electrodes all slow waves ( inlcudes NAN rows)
OnsetT.Properties.VariableNames=saveWaves.ChanLabels;
DurT.Properties.VariableNames=saveWaves.ChanLabels;
AmpT.Properties.VariableNames=saveWaves.ChanLabels;

FrontalElectrodes = ismember(saveWaves.ChanLabels,{'FPz'});
CentralElectrodes = ismember(saveWaves.ChanLabels,{'Cz'});
BackElectrodes = ismember(saveWaves.ChanLabels,{'Oz'});

% VERSION A: inlcude all waves
% % we need to get timecourse front
vOnset_F = reshape(table2array(OnsetT(:,FrontalElectrodes)), 1, [])';
vOnset_F(isnan(vOnset_F))=[];
vDur_F = reshape(table2array(DurT(:,FrontalElectrodes)), 1, [])';
vDur_F(isnan(vDur_F))=[];
vAmp_F= reshape(table2array(AmpT(:,FrontalElectrodes)), 1, [])';
vAmp_F(isnan(vAmp_F))=[];
temp_F = [vOnset_F, vDur_F, vAmp_F];
temp_F = sortrows(temp_F,[1,3],'ascend');
%center
vOnset_C = reshape(table2array(OnsetT(:,CentralElectrodes)), 1, [])';
vOnset_C(isnan(vOnset_C))=[];
vDur_C = reshape(table2array(DurT(:,CentralElectrodes)), 1, [])';
vDur_C(isnan(vDur_C))=[];
vAmp_C= reshape(table2array(AmpT(:,CentralElectrodes)), 1, [])';
vAmp_C(isnan(vAmp_C))=[];
temp_C = [vOnset_C, vDur_C, vAmp_C];
temp_C = sortrows(temp_C,[1,3],'ascend');
%back
vOnset_B = reshape(table2array(OnsetT(:,BackElectrodes)), 1, [])';
vOnset_B(isnan(vOnset_B))=[];
vDur_B = reshape(table2array(DurT(:,BackElectrodes)), 1, [])';
vDur_B(isnan(vDur_B))=[];
vAmp_B= reshape(table2array(AmpT(:,BackElectrodes)), 1, [])';
vAmp_B(isnan(vAmp_B))=[];
temp_B = [vOnset_B, vDur_B, vAmp_B];
temp_B = sortrows(temp_B,[1,3],'ascend'); % here we sort rows into ascending order based on onset and secondly based on amplitude

% Need better way but for now removing duplicate onsets (smaller amplitudes
% are removed
dx1 = find(diff(temp_F(:,1))==0);
temp_F(dx1,:)=[];
dx2 = find(diff(temp_C(:,1))==0);
temp_C(dx2,:)=[];
dx3 = find(diff(temp_B(:,1))==0);
temp_B(dx3,:)=[];

% recacluate into seconds
temp_F(:,1)=temp_F(:,1)/EEG.srate;
temp_F(:,2)=temp_F(:,2)/EEG.srate;
temp_C(:,1)=temp_C(:,1)/EEG.srate;
temp_C(:,2)=temp_C(:,2)/EEG.srate;
temp_B(:,1)=temp_B(:,1)/EEG.srate;
temp_B(:,2)=temp_B(:,2)/EEG.srate;

%into table
FTable = array2table(temp_F, 'VariableNames', {'Onset', 'Duration', 'Amplitude'});
CTable = array2table(temp_C, 'VariableNames', {'Onset', 'Duration', 'Amplitude'});
BTable = array2table(temp_B, 'VariableNames', {'Onset', 'Duration', 'Amplitude'});

%%
% SAVE 
savefilenameB1SF = ['swTimes_', subFolder '_rest_FPz.csv'];
writetable(FTable,[[TimingPath subFolder] filesep savefilenameB1SF]);
savefilenameB1SC = ['swTimes_', subFolder '_rest_Cz.csv'];
writetable(CTable,[[TimingPath subFolder] filesep savefilenameB1SC]);
savefilenameB1SB = ['swTimes_', subFolder '_rest_Oz.csv'];
writetable(BTable,[[TimingPath subFolder] filesep savefilenameB1SB]);

end