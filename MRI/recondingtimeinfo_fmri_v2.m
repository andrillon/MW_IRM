% Recoding trial and scanner trigger information into trial onset and
% duration for fMRI analysis


% The code creates timing files for further fMRI analysis. We need the
% behavioural data log files to get trial types (go or nogo) and probe
% responses. We also need the EEG marker file to get accurate information
% on the block start and first volume. 

%%
% INFO
% test_res
% 1: nblock 
% 2: this_blockcond 
% 3: thiset 
% 4: ntrial 
% 5: this_seq_trial
% 6: TargetID 
% 7: thisresp 
% 8: stimonset 
% 9: FlipSec-stim disappear
% 10: dur_face_presentation+dur_face_blank 
% 11: thisresptime(1)  
% 12:this_nogo 
% 13: this_go
%
% probe_res
% 1: this_probe 
% 2: this_probetime -- randomised probe time
% 3: startProbe --stop flip (next flip to randomised probe time)
% 4: endProbe -- after waitforscanner
% 5: nblock 
% 6: this_blockcond 
% 7: ntrial 
% 8-10: probe_responses -button
% 11-13: probe resp time
% 14-16: probe Q onset time
% 17-19: probe resp
%
% TR: 1560 ms
% TE: 30 ms
% Number slices: 42
% Number of volumes: depends on run
%
% Tiggers
% A: beginning of task
% E: end of task
% B: start of block (expected 4)
% C: end of block (expected 4)
% P: start of a probe
% Q: end of a probe
% S: MRI triggers
% T: start of a trial

% subject s205, s207, s208, s210, s214, s218, s219, s220, s224, s228, s239 needed to special selection for
% correct timings

%%
clear; 
close all;

path = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/';
behavFiles = dir([path 'Shared drives/MW_fMRI_EEG/Data/Behav/*.mat']);
behavFiles = behavFiles(2:end);

markerFiles = dir([path 'Shared drives/MW_fMRI_EEG/Data/EEG/Raw/*.vmrk']);
markerFiles = markerFiles([1:11,13:2:19,20:2:end],:);

for j = 3%9:size(behavFiles,1)

C = strsplit(behavFiles(j).name, '_');
fileToRead = [behavFiles(j).folder filesep  behavFiles(j).name];
TaskLog = load('-mat', fileToRead);

EEGMarkerFile = readlines([markerFiles(j).folder filesep  markerFiles(j).name]);

EEGMarkerFile = EEGMarkerFile(13:end-1,:); % remove info lines

%%
expectedDelimiters = 4;  % Number of expected delimiters
% Function to check delimiter count
validData = cellfun(@(x) sum(count(x, ',') == expectedDelimiters), EEGMarkerFile);

% Filter data to only include valid entries
EEGMarkerFile = EEGMarkerFile(validData==1,:);
Markers = split(EEGMarkerFile,","); 

Events = table;
Events.type = Markers(:,2);
Events.latency = str2double(Markers(:,3));

switch true
    case (strcmp(C{3}, 's205')==1)
        Events = Events(41:end,:); % we have some extra triggers to remove
    case (strcmp(C{3}, 's213')==1)
        Events = Events(60:end,:);% we have some extra triggers to remove
    case (strcmp(C{3}, 's214')==1)
        Events = Events(80:end,:);% we have some extra triggers to remove
    case (strcmp(C{3}, 's210')==1) % matlab crash, needed to restart task 2X so unrecorded sections are removed
        Events = Events([1:1500,2208:end],:);
    case (strcmp(C{3}, 's219')==1) % we have some extra triggers to remove
        Events = Events(41:end,:);
    case (strcmp(C{3}, 's228')==1)
        Events = Events(41:end,:);% we have some extra triggers to remove
    case (strcmp(C{3}, 's239')==1)
        Events = Events(41:end,:);% we have some extra triggers to remove
end


%%

% find rows in Events with give triggers
A1_list = [];
E1_list = [];
B1_list = [];
C1_list = [];
P1_list = [];
Q1_list = [];
T1_list = [];
S1_list = [];

k = 1;
for i = 1:size(Events,1)
    if strcmp (Events.type(i), 'B  1')
        B1_list(1,k) = Events.latency(i);
        k = k + 1;
    elseif strcmp (Events.type(i), 'A  1')
        A1_list(1,k) = Events.latency(i);
        k = k + 1;
    elseif strcmp (Events.type(i), 'C  1')
        C1_list(1,k) = Events.latency(i);
        k = k + 1;
    elseif strcmp (Events.type(i), 'E  1')
        E1_list(1,k) = Events.latency(i);
        k = k + 1;
    elseif strcmp (Events.type(i), 'S  1')
        S1_list(1,k) = Events.latency(i);
        k = k + 1; 
    elseif strcmp (Events.type(i), 'T  1')
        T1_list(1,k) = Events.latency(i);
        k = k + 1;
    elseif strcmp (Events.type(i), 'P  1')
        P1_list(1,k) = Events.latency(i);
        k = k + 1;
    elseif strcmp (Events.type(i), 'Q  1')
        Q1_list(1,k) = Events.latency(i);
        k = k + 1;
    end
    EEGTRIG.B1_list = nonzeros(B1_list);
    EEGTRIG.A1_list = nonzeros(A1_list);
    EEGTRIG.C1_list = nonzeros(C1_list);
    EEGTRIG.E1_list = nonzeros(E1_list);
    EEGTRIG.S1_list = nonzeros(S1_list);
    EEGTRIG.T1_list = nonzeros(T1_list);
    EEGTRIG.P1_list = nonzeros(P1_list);
    EEGTRIG.Q1_list = nonzeros(Q1_list);
end


%%
% Separating trials into 4 blocks 
firstBlockStart = EEGTRIG.B1_list(1);
secBlockStart = EEGTRIG.B1_list(2);
thirdBlockStart = EEGTRIG.B1_list(3);
forthBlockStart = EEGTRIG.B1_list(4);
forthBlockEnd = EEGTRIG.C1_list(4);

% elseif strcmp (C{3}, 's239')
%     secBlockStart = 7509401;
% end


ind1 = find(firstBlockStart<EEGTRIG.T1_list & EEGTRIG.T1_list<secBlockStart);
Block1_onsettimes = EEGTRIG.T1_list(ind1,:);
Block1_onsettimes = Block1_onsettimes(1:end-1,:); % because last trial trigger in the end of block in EEG data is invalid

ind2 = find(secBlockStart<EEGTRIG.T1_list & EEGTRIG.T1_list<thirdBlockStart);
Block2_onsettimes = EEGTRIG.T1_list(ind2,:);
Block2_onsettimes = Block2_onsettimes(1:end-1,:); 

ind3 = find(thirdBlockStart<EEGTRIG.T1_list & EEGTRIG.T1_list<forthBlockStart);
Block3_onsettimes = EEGTRIG.T1_list(ind3,:);
Block3_onsettimes = Block3_onsettimes(1:end-1,:); 

ind4 = find(forthBlockStart<EEGTRIG.T1_list & EEGTRIG.T1_list<forthBlockEnd);
Block4_onsettimes = EEGTRIG.T1_list(ind4,:);
Block4_onsettimes = Block4_onsettimes(1:end-1,:); 

%%

% Sometimes the scanner starts before B1 trigger and the first S1 tigger
% is before B1 so the next section 1) finds the first S1 trigger after
% the Block start (B1) and then checks if the preceeding S1 trigger is
% exactly 1 TR (1.56 sec = 7800 datapoints) before the S1 trigger -->
% if that's true the preceeding S1 is used as the first S1 trigger in the
% block.

S1ind1 = find(firstBlockStart<EEGTRIG.S1_list,1, "first");
S1ind2 = find(secBlockStart<EEGTRIG.S1_list,1, "first");
S1ind3 = find(thirdBlockStart<EEGTRIG.S1_list,1, "first");
S1ind4 = find(forthBlockStart<EEGTRIG.S1_list,1, "first");


s = 0;
while s == 0
    Block1_firstS1 = EEGTRIG.S1_list(S1ind1,:);
    if S1ind1<2
        s = 1;
    elseif S1ind1>=2
        if round(Block1_firstS1-EEGTRIG.S1_list(S1ind1-1,:),-2)==7800 %adding round to correct for potential 1 datapoint difference
        S1ind1 = S1ind1-1;
        else
        s = 1;
        end
    end
end



s = 0;
while s == 0
Block2_firstS1 = EEGTRIG.S1_list(S1ind2,:);
    if round(Block2_firstS1-EEGTRIG.S1_list(S1ind2-1,:),-2)==7800
        S1ind2 = S1ind2-1;
        if S1ind2<2
            Block2_firstS1 = EEGTRIG.S1_list(S1ind2,:);
            s=1;
        end 
    else
        s=1;
    end
end

s = 0;
while s == 0
Block3_firstS1 = EEGTRIG.S1_list(S1ind3,:);
    if round(Block3_firstS1-EEGTRIG.S1_list(S1ind3-1,:),-2)==7800
        S1ind3 = S1ind3-1;
    else
        s=1;
    end
end

s = 0;
while s == 0
Block4_firstS1 = EEGTRIG.S1_list(S1ind4,:);
    if round(Block4_firstS1-EEGTRIG.S1_list(S1ind4-1,:),-2)==7800
        S1ind4 = S1ind4-1;
    else
        s=1;
    end
end

%% Check difference between EEG triggers and task logs
% 
% TasklogTimeDiff1 = min(TaskLog.test_res(TaskLog.test_res(:,1)==1,8))-TaskLog.all_tstartblock(1,1);
% MarkerTimeDiff1 = ((min(Block1_onsettimes))-firstBlockStart)/5000;
% fprintf('Ppt: %s Block1\nTask time difference: %4.2f\nEEG marker time difference: %4.2f\n',C{3},TasklogTimeDiff1,MarkerTimeDiff1);
% 
% TasklogTimeDiff2 = min(TaskLog.test_res(TaskLog.test_res(:,1)==2,8))-TaskLog.all_tstartblock(1,2);
% MarkerTimeDiff2 = ((min(Block2_onsettimes))-secBlockStart)/5000;
% fprintf('Ppt: %s Block2\nTask time difference: %4.2f\nEEG marker time difference: %4.2f\n',C{3},TasklogTimeDiff2,MarkerTimeDiff2);
% 
% TasklogTimeDiff3 = min(TaskLog.test_res(TaskLog.test_res(:,1)==3,8))-TaskLog.all_tstartblock(1,3);
% MarkerTimeDiff3 = ((min(Block3_onsettimes))-thirdBlockStart)/5000;
% fprintf('Ppt: %s Block3\nTask time difference: %4.2f\nEEG marker time difference: %4.2f\n',C{3},TasklogTimeDiff3,MarkerTimeDiff3);
% 
% TasklogTimeDiff4 = min(TaskLog.test_res(TaskLog.test_res(:,1)==4,8))-TaskLog.all_tstartblock(1,4);
% MarkerTimeDiff4 = ((min(Block4_onsettimes))-forthBlockStart)/5000;
% fprintf('Ppt: %s Block4\nTask time difference: %4.2f\nEEG marker time difference: %4.2f\n',C{3},TasklogTimeDiff4,MarkerTimeDiff4);
% 

%%

% Subtract first S1 trigger time from trial onset times and transform them
% into milisecond

   % trial
        Block1_trial=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2_trial=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3_trial=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4_trial=(Block4_onsettimes-Block4_firstS1)/5000;

        block_trial = {Block1_trial, Block2_trial, Block3_trial, Block4_trial};

        for bln=1:length(block_trial)

            block_trial{bln}(:,2)=0.5;
            block_trial{bln}=array2table(block_trial{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_trial{bln},1)
                Trials = TaskLog.test_res(TaskLog.test_res(:,1)==bln,:);
                if Trials(tn,13)==1
                    block_trial{bln}.Type(tn) = {'ValidGo'};
                elseif Trials(tn,13)==0
                    block_trial{bln}.Type(tn) = {'Miss'};
                elseif isnan(Trials(tn,13)) && Trials(tn,12)==1
                    block_trial{bln}.Type(tn) = {'ValidNoGo'};
                elseif isnan(Trials(tn,13)) && Trials(tn,12)==0
                    block_trial{bln}.Type(tn) = {'FalseAlarm'};
                end
            end

        end

        % Response
        Block1_trialresp=((Block1_onsettimes-Block1_firstS1)/5000)+(TaskLog.test_res(TaskLog.test_res(:,1)==1,11)-TaskLog.test_res(TaskLog.test_res(:,1)==1,8));
        Block2_trialresp=((Block2_onsettimes-Block2_firstS1)/5000)+(TaskLog.test_res(TaskLog.test_res(:,1)==2,11)-TaskLog.test_res(TaskLog.test_res(:,1)==2,8));
        Block3_trialresp=((Block3_onsettimes-Block3_firstS1)/5000)+(TaskLog.test_res(TaskLog.test_res(:,1)==3,11)-TaskLog.test_res(TaskLog.test_res(:,1)==3,8));
        Block4_trialresp=((Block4_onsettimes-Block4_firstS1)/5000)+(TaskLog.test_res(TaskLog.test_res(:,1)==4,11)-TaskLog.test_res(TaskLog.test_res(:,1)==4,8));
        
        block_trialresp = {Block1_trialresp, Block2_trialresp, Block3_trialresp, Block4_trialresp};

        for bln=1:length(block_trialresp)

            block_trialresp{bln}(:,2)=0;
            block_trialresp{bln}=array2table(block_trialresp{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_trialresp{bln},1)
                Trials = TaskLog.test_res(TaskLog.test_res(:,1)==bln,:);
                if Trials(tn,13)==1
                    block_trialresp{bln}.Type(tn) = {'ValidGoResp'};
                elseif Trials(tn,13)==0
                    block_trialresp{bln}.Type(tn) = {'MissResp'}; %these should not exist
                elseif isnan(Trials(tn,13)) && Trials(tn,12)==1
                    block_trialresp{bln}.Type(tn) = {'ValidNoGoResp'}; %these should not exist
                elseif isnan(Trials(tn,13)) && Trials(tn,12)==0
                    block_trialresp{bln}.Type(tn) = {'FalseAlarmResp'};  
                end
            end

            block_trialresp{bln}(isnan(block_trialresp{bln}.Onset),:)=[];
            block_trialresp{bln}(1,:)=[];

        end

        %probes
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block2Diff=(secBlockStart-Block2_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(forthBlockStart-Block4_firstS1)/5000;
        
        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block2ScanStart = TaskLog.all_tstartblock(1,2)-Block2Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;

        % probe stop
        Block1_probestop=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,3)-Block1ScanStart;
        Block2_probestop=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,3)-Block2ScanStart;
        Block3_probestop=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,3)-Block3ScanStart;
        Block4_probestop=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,3)-Block4ScanStart;

        block_probestop = {Block1_probestop, Block2_probestop, Block3_probestop, Block4_probestop};

        for bln=1:length(block_probestop)

            block_probestop{bln}(:,2)=0.5;
            block_probestop{bln}=array2table(block_probestop{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probestop{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,17)==1
                    block_probestop{bln}.Type(tn) = {'ProbeStopONTASK'};
                elseif Trials(tn,17)==2
                    block_probestop{bln}.Type(tn) = {'ProbeStopMINDWANDER'};
                elseif Trials(tn,17)==3
                    block_probestop{bln}.Type(tn) = {'ProbeStopMINDBLANK'};
                elseif Trials(tn,17)==4
                    block_probestop{bln}.Type(tn) = {'ProbeStopNORECALL'};
                end
            end
        end

        % probe Q1
        Block1_probeq1=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,14)-Block1ScanStart;
        Block2_probeq1=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,14)-Block2ScanStart;
        Block3_probeq1=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,14)-Block3ScanStart;
        Block4_probeq1=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,14)-Block4ScanStart;

        block_probeq1 = {Block1_probeq1, Block2_probeq1, Block3_probeq1, Block4_probeq1};

        for bln=1:length(block_probeq1)
 
            block_probeq1{bln}(:,2)=TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,11)-TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,14);
            block_probeq1{bln}=array2table(block_probeq1{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeq1{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,17)==1
                    block_probeq1{bln}.Type(tn) = {'Q1ONTASK'};
                elseif Trials(tn,17)==2
                    block_probeq1{bln}.Type(tn) = {'Q1MINDWANDER'};
                elseif Trials(tn,17)==3
                    block_probeq1{bln}.Type(tn) = {'Q1MINDBLANK'};
                elseif Trials(tn,17)==4
                    block_probeq1{bln}.Type(tn) = {'Q1NORECALL'};
                end
            end
        end

        % probe Q2
        Block1_probeq2=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,15)-Block1ScanStart;
        Block2_probeq2=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,15)-Block2ScanStart;
        Block3_probeq2=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,15)-Block3ScanStart;
        Block4_probeq2=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,15)-Block4ScanStart;

        block_probeq2 = {Block1_probeq2, Block2_probeq2, Block3_probeq2, Block4_probeq2};

        for bln=1:length(block_probeq2)
 
            block_probeq2{bln}(:,2)=TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,12)-TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,15);
            block_probeq2{bln}=array2table(block_probeq2{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeq2{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,18)==1
                    block_probeq2{bln}.Type(tn) = {'Q2MWroom'};
                elseif Trials(tn,18)==2
                    block_probeq2{bln}.Type(tn) = {'Q2MWpersonal'};
                elseif Trials(tn,18)==3
                    block_probeq2{bln}.Type(tn) = {'Q2MWtask'};
                else
                    block_probeq2{bln}.Type(tn) = {'NaN'};
                end
            end
            block_probeq2{bln}(isnan(block_probeq2{bln}.Onset),:)=[];
        end

        % probe Q3
        Block1_probeq3=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,16)-Block1ScanStart;
        Block2_probeq3=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,16)-Block2ScanStart;
        Block3_probeq3=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,16)-Block3ScanStart;
        Block4_probeq3=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,16)-Block4ScanStart;

        block_probeq3 = {Block1_probeq3, Block2_probeq3, Block3_probeq3, Block4_probeq3};

        for bln=1:length(block_probeq3)
 
            block_probeq3{bln}(:,2)=TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,13)-TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,16);
            block_probeq3{bln}=array2table(block_probeq3{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeq3{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,19)==1
                    block_probeq3{bln}.Type(tn) = {'Q3ExAlert'};
                elseif Trials(tn,19)==2
                    block_probeq3{bln}.Type(tn) = {'Q3Alert'};
                elseif Trials(tn,19)==3
                    block_probeq3{bln}.Type(tn) = {'Q3Sleepy'};
                elseif Trials(tn,19)==4
                    block_probeq3{bln}.Type(tn) = {'Q3ExSleepy'};
                end
            end
        end

        % probe Q1resp
        Block1_probeq1resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,11)-Block1ScanStart;
        Block2_probeq1resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,11)-Block2ScanStart;
        Block3_probeq1resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,11)-Block3ScanStart;
        Block4_probeq1resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,11)-Block4ScanStart;

        block_probeq1resp = {Block1_probeq1resp, Block2_probeq1resp, Block3_probeq1resp, Block4_probeq1resp};

        for bln=1:length(block_probeq1resp)
 
            block_probeq1resp{bln}(:,2)=0;
            block_probeq1resp{bln}=array2table(block_probeq1resp{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeq1resp{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,17)==1
                    block_probeq1resp{bln}.Type(tn) = {'Q1RespONTASK'};
                elseif Trials(tn,17)==2
                    block_probeq1resp{bln}.Type(tn) = {'Q1RespMINDWANDER'};
                elseif Trials(tn,17)==3
                    block_probeq1resp{bln}.Type(tn) = {'Q1RespMINDBLANK'};
                elseif Trials(tn,17)==4
                    block_probeq1resp{bln}.Type(tn) = {'Q1RespNORECALL'};
                end
            end
        end

        % probe Q2resp
        Block1_probeq2resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,12)-Block1ScanStart;
        Block2_probeq2resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,12)-Block2ScanStart;
        Block3_probeq2resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,12)-Block3ScanStart;
        Block4_probeq2resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,12)-Block4ScanStart;

        block_probeq2resp = {Block1_probeq2resp, Block2_probeq2resp, Block3_probeq2resp, Block4_probeq2resp};

        for bln=1:length(block_probeq2resp)
 
            block_probeq2resp{bln}(:,2)=0;
            block_probeq2resp{bln}=array2table(block_probeq2resp{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeq2resp{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,18)==1
                    block_probeq2resp{bln}.Type(tn) = {'Q2RespMWroom'};
                elseif Trials(tn,18)==2
                    block_probeq2resp{bln}.Type(tn) = {'Q2RespMWpersonal'};
                elseif Trials(tn,18)==3
                    block_probeq2resp{bln}.Type(tn) = {'Q2RespMWtask'};
                else
                    block_probeq2resp{bln}.Type(tn) = {'NaN'};
                end
            end
            block_probeq2resp{bln}(isnan(block_probeq2resp{bln}.Onset),:)=[];
        end

        % probe Q3resp
        Block1_probeq3resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,13)-Block1ScanStart;
        Block2_probeq3resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,13)-Block2ScanStart;
        Block3_probeq3resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,13)-Block3ScanStart;
        Block4_probeq3resp=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,13)-Block4ScanStart;

        block_probeq3resp = {Block1_probeq3resp, Block2_probeq3resp, Block3_probeq3resp, Block4_probeq3resp};

        for bln=1:length(block_probeq3resp)
 
            block_probeq3resp{bln}(:,2)=0;
            block_probeq3resp{bln}=array2table(block_probeq3resp{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeq3resp{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,19)==1
                    block_probeq3resp{bln}.Type(tn) = {'Q3RespExAlert'};
                elseif Trials(tn,19)==2
                    block_probeq3resp{bln}.Type(tn) = {'Q3RespAlert'};
                elseif Trials(tn,19)==3
                    block_probeq3resp{bln}.Type(tn) = {'Q3RespSleepy'};
                elseif Trials(tn,19)==4
                    block_probeq3resp{bln}.Type(tn) = {'Q3RespExSleepy'};
                end
            end
        end

        % probe end
        Block1_probeend=TaskLog.probe_res(TaskLog.probe_res(:,5)==1,13)-Block1ScanStart;
        Block2_probeend=TaskLog.probe_res(TaskLog.probe_res(:,5)==2,13)-Block2ScanStart;
        Block3_probeend=TaskLog.probe_res(TaskLog.probe_res(:,5)==3,13)-Block3ScanStart;
        Block4_probeend=TaskLog.probe_res(TaskLog.probe_res(:,5)==4,13)-Block4ScanStart;

        block_probeend = {Block1_probeend, Block2_probeend, Block3_probeend, Block4_probeend};

        for bln=1:length(block_probeend)

            block_probeend{bln}(:,2)=TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,4)-TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,13);
            block_probeend{bln}=array2table(block_probeend{bln}, 'VariableNames',{'Onset','Duration'});

            for tn=1:size(block_probeend{bln},1)
                Trials = TaskLog.probe_res(TaskLog.probe_res(:,5)==bln,:);
                if Trials(tn,17)==1
                    block_probeend{bln}.Type(tn) = {'ProbeEndONTASK'};
                elseif Trials(tn,17)==2
                    block_probeend{bln}.Type(tn) = {'ProbeEndMINDWANDER'};
                elseif Trials(tn,17)==3
                    block_probeend{bln}.Type(tn) = {'ProbeEndMINDBLANK'};
                elseif Trials(tn,17)==4
                    block_probeend{bln}.Type(tn) = {'ProbeEndNORECALL'};
                end
            end
        end
    

        % Merge tables 
        Block1 = [block_trial{1}; block_trialresp{1}; block_probestop{1}; block_probeq1{1}; block_probeq1resp{1}; block_probeq2{1}; block_probeq2resp{1}; block_probeq3{1}; block_probeq3resp{1}; block_probeend{1}];
        Block2 = [block_trial{2}; block_trialresp{2}; block_probestop{2}; block_probeq1{2}; block_probeq1resp{2}; block_probeq2{2}; block_probeq2resp{2}; block_probeq3{2}; block_probeq3resp{2}; block_probeend{2}];
        Block3 = [block_trial{3}; block_trialresp{3}; block_probestop{3}; block_probeq1{3}; block_probeq1resp{3}; block_probeq2{3}; block_probeq2resp{3}; block_probeq3{3}; block_probeq3resp{3}; block_probeend{3}];
        Block4 = [block_trial{4}; block_trialresp{4}; block_probestop{4}; block_probeq1{4}; block_probeq1resp{4}; block_probeq2{4}; block_probeq2resp{4}; block_probeq3{4}; block_probeq3resp{4}; block_probeend{4}];

        % Sort the table by 'Onset'
        Block1 = sortrows(Block1, 'Onset');
        Block2 = sortrows(Block2, 'Onset');
        Block3 = sortrows(Block3, 'Onset');
        Block4 = sortrows(Block4, 'Onset');

if (strcmp(C{3}, 's207')==1) %we needed to stop the scanner in the middle of Block 2
    Block22 = Block2;
    Block2 = Block2(1:352,:);
    Block22 = Block22(353:end,:);
    Block22Diff=((5086057-7848712)/5000);
    Block22.Onset = Block22.Onset+Block22Diff;
    Block22 = Block22(Block22.Onset>0,:);
elseif (strcmp(C{3}, 's208')==1) %ppt did not take a break between block 1 and 2, some trials are missing from block 2
    Block1 = [Block1; Block2];
    scanEnd = (7766852-Block1_firstS1)/5000;
    Block1 = Block1(Block1.Onset<scanEnd,:); 
    secBlockStart2 = (8704388-Block1_firstS1)/5000;
    Block2 = Block2(Block2.Onset>secBlockStart2,:);
    Block2.Onset = Block2.Onset-secBlockStart2;
elseif (strcmp(C{3}, 's214')==1) %Block4 started before scanning started so some trials need to be removed
    Block4Diff=(14258045-Block4_firstS1)/5000;
    Block4.Onset = Block4.Onset-Block4Diff;
    Block4 = Block4(Block4.Onset>0,:);
elseif (strcmp(C{3}, 's218')==1) %ppt did not take a break between block 1 and 2
    Block1 = [Block1; Block2];
    scanEnd = (7742964-Block1_firstS1)/5000;
    Block1 = Block1(Block1.Onset<scanEnd,:); 
    secBlockStart2 = (8410022-Block1_firstS1)/5000;
    Block2 = Block2(Block2.Onset>secBlockStart2,:);
    Block2.Onset = Block2.Onset-secBlockStart2;
elseif (strcmp(C{3}, 's220')==1) %Block4 started before scanning started so some trials need to be removed
    Block4Diff=(13483947-Block4_firstS1)/5000;
    Block4.Onset = Block4.Onset-Block4Diff;
elseif (strcmp(C{3}, 's228')==1) %ppt did not take a break between block 3 and 4, some trials are missing from block 4
    Block3 = [Block3; Block4];
    scanEnd = (15965961-Block3_firstS1)/5000;
    Block3 = Block3(Block3.Onset<scanEnd,:); 
    forthBlockStart2 = (16890516-Block3_firstS1)/5000;
    Block4 = Block4(Block4.Onset>forthBlockStart2,:);
    Block4.Onset = Block4.Onset-forthBlockStart2;
elseif (strcmp(C{3}, 's239')==1) %Block2 started before scanning started
    Block2Diff=(7537490-Block2_firstS1)/5000;
    Block2.Onset = Block2.Onset-Block2Diff;

end


% Catch response errors
if any (Block1.Duration< 0)
   fprintf('Timing inconsitency in Task Block 1 for %s in %d instances\n',C{3}, sum(Block1.Duration< 0));
end

if any (Block2.Duration< 0)
   fprintf('Timing inconsitency in Task Block 2 for %s in %d instances\n',C{3}, sum(Block2.Duration< 0));
end 

if any (Block3.Duration< 0)
   fprintf('Timing inconsitency in Task Block 3 for %s in %d instances\n',C{3}, sum(Block3.Duration< 0));
end

if any (Block4.Duration< 0)
   fprintf('Timing inconsitency in Task Block 4 for %s in %d instances\n',C{3}, sum(Block4.Duration< 0));
end

% Catch negative onset times
if any (Block1.Onset< 0) || any (Block2.Onset< 0) || any (Block3.Onset< 0) || any (Block4.Onset< 0)
   %fprintf('Negative Onset times for %s\n\n',C{3}); %s224; s212
   %break
   Block1 = Block1(Block1.Onset>0,:);
   Block2 = Block2(Block2.Onset>0,:);
   Block3 = Block3(Block3.Onset>0,:);
   Block4 = Block4(Block4.Onset>0,:);
end

% % Catch late start
% if Block1.Onset(1,1)> 20 || Block2.Onset(1,1)> 20 || Block3.Onset(1,1)> 20 || Block4.Onset(1,1)> 20
%    fprintf('Non-zero start times for %s\n\n',C{3}); 
%    %break
% end




%% SAVE 

%savepath =[path 'Shared drives/MW_fMRI_EEG/Data/MRI/TimingforfMRI/'];
%savepath =[path 'Shared drives/MW_fMRI_EEG/Data/MRI/Timing2'];

% % create folder
% if exist([savepath C{3}]) == 0
%     fprintf('Folder does not exist. Creating: %s \n',C{3})
%     mkdir([savepath C{3}])
%     mkdir([savepath C{3} '/Block1'])
%     mkdir([savepath C{3} '/Block2'])
%     mkdir([savepath C{3} '/Block3'])
%     mkdir([savepath C{3} '/Block4'])
% else
%     fprintf('%s Folder exists. \n',C{3})
% end

savefilenameB1 = ['fulldesignmat_', C{3} '_Block1.csv'];
savefilenameB2 = ['fulldesignmat_', C{3} '_Block2.csv'];
savefilenameB3 = ['fulldesignmat_', C{3} '_Block3.csv'];
savefilenameB4 = ['fulldesignmat_', C{3} '_Block4.csv'];


% writetable(Block1,[[savepath C{3} '/Block1'] filesep savefilenameB1]);
% writetable(Block2,[[savepath C{3} '/Block2'] filesep savefilenameB2]);
% writetable(Block3,[[savepath C{3} '/Block3'] filesep savefilenameB3]);
% writetable(Block4,[[savepath C{3} '/Block4'] filesep savefilenameB4]);
% 
% 
% if  strcmp (C{3}, 's207')
     savefilenameB22 = ['fulldesignmat_', C{3} '_Block22.csv'];
% %     writetable(Block22,[[savepath C{3} '/Block2-2'] filesep savefilenameB22]);
     writetable(Block22, [savepath filesep savefilenameB22]);  
% end
% 
% 
% writetable(Block1,[savepath filesep savefilenameB1]);
% writetable(Block2,[savepath filesep savefilenameB2]);
% writetable(Block3,[savepath filesep savefilenameB3]);
% writetable(Block4,[savepath filesep savefilenameB4]);

end

