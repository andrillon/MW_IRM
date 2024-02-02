function [FTable, CTable, BTable] = SW_TimeRecodingforMR (BlockSW, BlockP, saveWaves) 

% % create a tabel on onsets,durations,amplitudes
OnsetT = table;
DurT = table;
AmpT = table;


nPlist=unique(BlockSW(:,16))';
for nP=nPlist % for every probe 
    

    allOnset ={};
    allDuration ={};
    allAmplitude ={};

        
    for nEl=unique(BlockSW(:,3))' %and every electrode 
        onset=[];
        duration =[];
        amplitude = [];
        temp_waves=BlockSW(BlockSW(:,16)==nP & BlockSW(:,3)==nEl,:);% take all slow waves
        for nW=1:size(temp_waves,1)
            onset = [onset; BlockP(BlockP(:,1)==nP,2) + (temp_waves(nW,5)/1000)];%subtract negative onset value from probe onset in seconds
            duration = [duration; (temp_waves(nW,8)-temp_waves(nW,5))/1000]; % calculate duration in seconds
            amplitude = [amplitude; abs(temp_waves(nW,9))]; % take absolute amplitude
        end
        allOnset{nEl}=onset; % for each electrode, we then compy all onset times 
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
end

% outcome matrix: 63 Electrodes all slow waves ( inlcudes NAN rows)
OnsetT.Properties.VariableNames=saveWaves.ChanLabels;
DurT.Properties.VariableNames=saveWaves.ChanLabels;
AmpT.Properties.VariableNames=saveWaves.ChanLabels;

% select subset of electrodes
FrontalElectrodes = ismember(saveWaves.ChanLabels,{'F1','Fz', 'F2', 'FC1', 'FC2'});
CentralElectrodes = ismember(saveWaves.ChanLabels,{'C1','Cz', 'C2', 'CP1', 'CPz', 'CP2'});
BackElectrodes = ismember(saveWaves.ChanLabels,{'P1','Pz', 'P2', 'PO3', 'POz', 'PO4'});

% we need to get timecourse 

vOnset_F = reshape(table2array(OnsetT(:,FrontalElectrodes)), 1, [])';
vOnset_F(isnan(vOnset_F))=[];
vDur_F = reshape(table2array(DurT(:,FrontalElectrodes)), 1, [])';
vDur_F(isnan(vDur_F))=[];
vAmp_F= reshape(table2array(AmpT(:,FrontalElectrodes)), 1, [])';
vAmp_F(isnan(vAmp_F))=[];

temp_F = [vOnset_F, vDur_F, vAmp_F];
temp_F = sortrows(temp_F,[1,3],'ascend');

vOnset_C = reshape(table2array(OnsetT(:,CentralElectrodes)), 1, [])';
vOnset_C(isnan(vOnset_C))=[];
vDur_C = reshape(table2array(DurT(:,CentralElectrodes)), 1, [])';
vDur_C(isnan(vDur_C))=[];
vAmp_C= reshape(table2array(AmpT(:,CentralElectrodes)), 1, [])';
vAmp_C(isnan(vAmp_C))=[];

temp_C = [vOnset_C, vDur_C, vAmp_C];
temp_C = sortrows(temp_C,[1,3],'ascend');


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

FTable = array2table(temp_F, 'VariableNames', {'Onset', 'Duration', 'Amplitude'});
CTable = array2table(temp_C, 'VariableNames', {'Onset', 'Duration', 'Amplitude'});
BTable = array2table(temp_B, 'VariableNames', {'Onset', 'Duration', 'Amplitude'});


end