function [Table] = SW_TimeRecodingforMR_singlechan(BlockSW, BlockP, saveWaves, nchan) 

% set parameters to extract
onsetp = 8; % 5=start neg. wave; 8=peak neg wave;6=start pos wave; 10=peak pos wave
durp = 6; %7=end of wave
ampp = 9; % 9=neg peak, 11=pos peak, 4=peak2peak


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

Electrode = ismember(saveWaves.ChanLabels,saveWaves.ChanLabels(nchan));

% VERSION A: inlcude all waves
% % we need to get timecourse 
vOnset_F = reshape(table2array(OnsetT(:,Electrode)), 1, [])';
vOnset_F(isnan(vOnset_F))=[];
vDur_F = reshape(table2array(DurT(:,Electrode)), 1, [])';
vDur_F(isnan(vDur_F))=[];
vAmp_F= reshape(table2array(AmpT(:,Electrode)), 1, [])';
vAmp_F(isnan(vAmp_F))=[];
% keep reshaping the probe table 
vSWP_F =[];
vSWP_F = repmat(table2array(SWProbeT), sum(Electrode),1); 
tempOnset_F = reshape(table2array(OnsetT(:,Electrode)), 1, [])';
vSWP_F(isnan(tempOnset_F),:)=[];

temp_ch = [vOnset_F, vDur_F, vAmp_F, vSWP_F];
temp_ch = sortrows(temp_ch,[1,3],'ascend');

Table = array2table(temp_ch, 'VariableNames', {'Onset', 'Duration', 'Amplitude', 'ProbeNum', 'Q1', 'Q2', 'Q3'});

end