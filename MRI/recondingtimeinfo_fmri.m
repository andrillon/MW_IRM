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
% 2:this_blockcond 
% 3: thiset 
% 4: ntrial 
% 5: this_seq_trial
% 6: TargetID 
% 7: thisresp 
% 8: stimonset 
% 9: FlipSec 
% 10: dur_face_presentation+dur_face_blank 
% 11: thisresptime(1)  
% 12:this_nogo 
% 13: this_go
%
% probe_res
% 1: this_probe 
% 2: this_probetime 
% 3: startProbe 
% 4: endProbe 
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

% subject s205, s207, s208, s210, s214, s219, s220, s228, s239 needed to special selection for
% correct timings

%s205 Synch Off!
%%
clear; 
close all;

path = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/';

%behavFiles = dir('/Volumes/TOSHIBA/MW_EEGFMRI_PROJECT/EEGFMRI_BehavData');
behavFiles = dir([path 'Shared drives/MW_fMRI_EEG/Data/Behav/*.mat']);

behavFiles = behavFiles(2:end);

%markerFiles = dir('/Volumes/TOSHIBA/MW_EEGFMRI_PROJECT/EEGFMRI_EEGData/*vmrk');
markerFiles = dir([path 'Shared drives/MW_fMRI_EEG/Data/EEG/Raw/*.vmrk']);

markerFiles = markerFiles([44:2:end],:);

%markerFiles = markerFiles([18,20,22,23,25,27,29,31,33,35,37,39,41,43],:);
%markerFiles = markerFiles([2,4,5,7,10],:);

for j = 5 %% 1:size(behavFiles,1)

C = strsplit(behavFiles(j).name, '_');
fileToRead = [behavFiles(j).folder filesep  behavFiles(j).name];
TaskLog = load('-mat', fileToRead);

EEGMarkerFile = readlines([markerFiles(j).folder filesep  markerFiles(j).name]);


EEGMarkerFile = EEGMarkerFile(13:end-1,:);

%%
% Catch errors

try  
    Markers = split(EEGMarkerFile,","); 
catch ME 
    EEGMarkerFile = EEGMarkerFile(56:end,:);
    try
        Markers = split(EEGMarkerFile,",");
    catch ME
        EEGMarkerFile = EEGMarkerFile(18:end,:);
        try
            Markers = split(EEGMarkerFile,",");
        catch ME
            %EEGMarkerFile = EEGMarkerFile(1:5554,:);
            EEGMarkerFile = EEGMarkerFile(1:end-1,:);
            Markers = split(EEGMarkerFile,",");
        end
    end
end


Events = table;
Events.type = Markers(:,2);
Events.latency = str2double(Markers(:,3));



switch true
    case (strcmp(C{3}, 's205')==1)
        Events = Events(80:end,:);
    case (strcmp(C{3}, 's214')==1)
        Events = Events(80:end,:);
    case (strcmp(C{3}, 's210')==1)
        Events = Events([1:1500,2208:end],:);
    case (strcmp(C{3}, 's219')==1)
        Events = Events(88:end,:);
    case (strcmp(C{3}, 's228')==1)
        Events = Events(88:end,:);
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
% Catch ERROR related to B1 triggers

if length(EEGTRIG.B1_list) ~= 4 
    fprintf('Number of blocks in EEG data in %s do not add up to 4. Execution halted \n',C{3});
    continue;
else
    fprintf('Number of blocks in EEG data in %s is 4. Continue. \n',C{3});
end 

% We find the first S1 trigger after each block start and that will be
% zero, with that we can set up the time difference between T1 trials and 

%%
% Separating trials into 4 blocks 
firstBlockStart = EEGTRIG.B1_list(1);
secBlockStart = EEGTRIG.B1_list(2);
thirdBlockStart = EEGTRIG.B1_list(3);
forthBlockStart = EEGTRIG.B1_list(4);
forthBlockEnd = EEGTRIG.C1_list(4);

if  strcmp (C{3}, 's208')
    secBlockStart = 8704388;
elseif strcmp (C{3}, 's218')
    secBlockStart = 8415401;
elseif strcmp (C{3}, 's220')
    forthBlockStart = 13471201;
elseif strcmp (C{3}, 's228')
    forthBlockStart = 16885201;
elseif strcmp (C{3}, 's239')
    secBlockStart = 7509401;
end


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


if  strcmp (C{3}, 's207')
    Block22_onsettimes = Block2_onsettimes;
    Block2_onsettimes = Block2_onsettimes(1:176,:); 
    Block22_onsettimes = Block22_onsettimes(177:end,:); 
elseif strcmp (C{3}, 's228')
    ind3 = find(thirdBlockStart<EEGTRIG.T1_list & EEGTRIG.T1_list<15965961);
    Block3_onsettimes = EEGTRIG.T1_list(ind3,:);
    Block3_onsettimes = Block3_onsettimes(1:end-1,:); 

end

%%

% Sometimes the scanner starts before B1 trigger and the first S1 tigger
% is before B1 so the next section 1) finds the first S1 trigger after
% the Block start (B1) and then checks if the preceeding S1 trigger is
% exactly 1 TR (1.56 sec ~ 7800 datapoints) before the S1 trigger -->
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
        if Block1_firstS1-EEGTRIG.S1_list(S1ind1-1,:)==7800
        S1ind1 = S1ind1-1;
        else
        s = 1;
        end
    end
end



s = 0;
while s == 0
Block2_firstS1 = EEGTRIG.S1_list(S1ind2,:);
    if Block2_firstS1-EEGTRIG.S1_list(S1ind2-1,:)==7800
        S1ind2 = S1ind2-1;
     % Catch ERROR related to ind < 2
        if S1ind2<2
            fprintf('Index number in %s is invalid. Execution halted \n',C{3});
            s=1;
            continue;
        end 
    else
        s=1;
    end
end

s = 0;
while s == 0
Block3_firstS1 = EEGTRIG.S1_list(S1ind3,:);
    if Block3_firstS1-EEGTRIG.S1_list(S1ind3-1,:)==7800
        S1ind3 = S1ind3-1;
    else
        s=1;
    end
end

s = 0;
while s == 0
Block4_firstS1 = EEGTRIG.S1_list(S1ind4,:);
    if Block4_firstS1-EEGTRIG.S1_list(S1ind4-1,:)==7800
        S1ind4 = S1ind4-1;
    else
        s=1;
    end
end

%%

% Deduct first S1 trigger time from trial onset times and transform them
% into milisecond

switch true
    case (strcmp(C{3}, 's207')==1)
        Block1=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4=(Block4_onsettimes-Block4_firstS1)/5000;

        Block1 = [Block1,TaskLog.test_res(TaskLog.test_res(:,1)==1, [4,10,11,12,13])];
        Block1(:,4) = Block1(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==1,8);
        Block3 = [Block3,TaskLog.test_res(TaskLog.test_res(:,1)==3, [4,10,11,12,13])];
        Block3(:,4) = Block3(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==3,8);
        Block4 = [Block4,TaskLog.test_res(TaskLog.test_res(:,1)==4, [4,10,11,12,13])];
        Block4(:,4) = Block4(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==4,8);

        Block22=(Block22_onsettimes-7848712)/5000;
        fullBlock = TaskLog.test_res(TaskLog.test_res(:,1)==2, [4,10,11,12,13]);
        fullBlock(:,3) =fullBlock(:,3)-TaskLog.test_res(TaskLog.test_res(:,1)==2,8);
        Block2 = [Block2,fullBlock(1:176,:)];
        Block22 = [Block22,fullBlock(177:end,:)];
        Block22 = array2table(Block22,...
        'VariableNames',{'TrialOnsetTime','TrialNum','StimDur', 'RespRT', 'ThisNogo', 'ThisGo'});
    case (strcmp(C{3}, 's208')==1)
        Block1=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4=(Block4_onsettimes-Block4_firstS1)/5000;

        fullBlock = TaskLog.test_res(TaskLog.test_res(:,1)<3, [4,10,11,12,13]);
        fullBlock(:,3) =fullBlock(:,3)-TaskLog.test_res(TaskLog.test_res(:,1)<3,8);
        Block1 = [Block1,fullBlock(1:size(Block1,1),:)];
        Block2 = [Block2,fullBlock(size(Block1,1)+1:end,:)];
        scanEnd = (7766852-Block1_firstS1)/5000;
        Block1 = Block1(Block1(:,1)<scanEnd,:); 
        Block3 = [Block3,TaskLog.test_res(TaskLog.test_res(:,1)==3, [4,10,11,12,13])];
        Block3(:,4) = Block3(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==3,8);
        Block4 = [Block4,TaskLog.test_res(TaskLog.test_res(:,1)==4, [4,10,11,12,13])];
        Block4(:,4) = Block4(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==4,8);
    case (strcmp(C{3}, 's214')==1)
        Block1=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4=(Block4_onsettimes-14258045)/5000;

        Block1 = [Block1,TaskLog.test_res(TaskLog.test_res(:,1)==1, [4,10,11,12,13])];
        Block1(:,4) = Block1(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==1,8);
        Block2 = [Block2,TaskLog.test_res(TaskLog.test_res(:,1)==2, [4,10,11,12,13])];
        Block2(:,4) = Block2(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==2,8);
        Block3 = [Block3,TaskLog.test_res(TaskLog.test_res(:,1)==3, [4,10,11,12,13])];
        Block3(:,4) = Block3(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==3,8);
        Block4 = [Block4,TaskLog.test_res(TaskLog.test_res(:,1)==4, [4,10,11,12,13])];
        Block4(:,4) = Block4(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==4,8);
        Block4=Block4(25:end,:);

    case (strcmp(C{3}, 's218')==1)
        Block1=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4=(Block4_onsettimes-Block4_firstS1)/5000;

        fullBlock = TaskLog.test_res(TaskLog.test_res(:,1)<3, [4,10,11,12,13]);
        fullBlock(:,3) =fullBlock(:,3)-TaskLog.test_res(TaskLog.test_res(:,1)<3,8);

        Block1 = [Block1,fullBlock(1:size(Block1,1),:)];
        Block2 = [Block2,fullBlock(size(Block1,1)+1:end,:)];
        scanEnd = (7750764-Block1_firstS1)/5000;
        Block1 = Block1(Block1(:,1)<scanEnd,:); 

        Block3 = [Block3,TaskLog.test_res(TaskLog.test_res(:,1)==3, [4,10,11,12,13])];
        Block3(:,4) = Block3(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==3,8);
        Block4 = [Block4,TaskLog.test_res(TaskLog.test_res(:,1)==4, [4,10,11,12,13])];
        Block4(:,4) = Block4(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==4,8);

    case (strcmp(C{3}, 's228')==1)
        Block1=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4=(Block4_onsettimes-Block4_firstS1)/5000;

        Block1 = [Block1,TaskLog.test_res(TaskLog.test_res(:,1)==1, [4,10,11,12,13])];
        Block1(:,4) = Block1(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==1,8);
        Block2 = [Block2,TaskLog.test_res(TaskLog.test_res(:,1)==2, [4,10,11,12,13])];
        Block2(:,4) = Block2(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==2,8);

        fullBlock = TaskLog.test_res(TaskLog.test_res(:,1)>=3, [4,10,11,12,13]);
        fullBlock(:,3) =fullBlock(:,3)-TaskLog.test_res(TaskLog.test_res(:,1)>=3,8);

        Block3 = [Block3,fullBlock(1:size(Block3,1),:)];
        scanEnd = (15965961-Block3_firstS1)/5000;
        Block3 = Block3(Block3(:,1)<scanEnd,:); 
        Block4 = [Block4,fullBlock(end-size(Block4,1)+1:end,:)];

        
    otherwise
        Block1=(Block1_onsettimes-Block1_firstS1)/5000;
        Block2=(Block2_onsettimes-Block2_firstS1)/5000;
        Block3=(Block3_onsettimes-Block3_firstS1)/5000;
        Block4=(Block4_onsettimes-Block4_firstS1)/5000;

        Block1 = [Block1,TaskLog.test_res(TaskLog.test_res(:,1)==1, [4,10,11,12,13])];
        Block1(:,4) = Block1(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==1,8);
        Block2 = [Block2,TaskLog.test_res(TaskLog.test_res(:,1)==2, [4,10,11,12,13])];
        Block2(:,4) = Block2(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==2,8);
        Block3 = [Block3,TaskLog.test_res(TaskLog.test_res(:,1)==3, [4,10,11,12,13])];
        Block3(:,4) = Block3(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==3,8);
        Block4 = [Block4,TaskLog.test_res(TaskLog.test_res(:,1)==4, [4,10,11,12,13])];
        Block4(:,4) = Block4(:,4)-TaskLog.test_res(TaskLog.test_res(:,1)==4,8);
end

% Catch late scanning start
if any (Block1(:,1)< 0)
   fprintf('Timing inconsitency in Task Block 1 for %s \n',C{3});
   Block1 = Block1(Block1(:,1)>0, :);
end

if any (Block2(:,1)< 0)
   fprintf('Timing inconsitency in Task Block 2 for %s \n',C{3});
   Block2 = Block2(Block2(:,1)>0, :);
end 

if any (Block3(:,1)< 0)
   fprintf('Timing inconsitency in Task Block 3 for %s \n',C{3});
   Block3 = Block3(Block3(:,1)>0, :);
end

if any (Block4(:,1)< 0)
   fprintf('Timing inconsitency in Task Block 4 for %s \n',C{3});
   Block4 = Block4(Block4(:,1)>0, :);
end

Block1 = array2table(Block1,...
    'VariableNames',{'TrialOnsetTime','TrialNum','StimDur', 'RespRT', 'ThisNogo', 'ThisGo'});
Block2 = array2table(Block2,...
    'VariableNames',{'TrialOnsetTime','TrialNum','StimDur', 'RespRT', 'ThisNogo', 'ThisGo'});
Block3 = array2table(Block3,...
    'VariableNames',{'TrialOnsetTime','TrialNum','StimDur', 'RespRT', 'ThisNogo', 'ThisGo'});
Block4 = array2table(Block4,...
    'VariableNames',{'TrialOnsetTime','TrialNum','StimDur', 'RespRT', 'ThisNogo', 'ThisGo'});


%% Probes

switch true
    case (strcmp(C{3}, 's207')==1)
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block2Diff=(secBlockStart-Block2_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(forthBlockStart-Block4_firstS1)/5000;
        Block22Diff=(secBlockStart-7848712)/5000;

        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block2ScanStart = TaskLog.all_tstartblock(1,2)-Block2Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;
        Block22ScanStart = TaskLog.all_tstartblock(1,2)-Block22Diff;

        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==1, [1,3,4,7,17:19, 16]);
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block22Probe = Block2Probe;
        Block2Probe = Block2Probe(Block2Probe(:,2)<Block22ScanStart,:);
        Block22Probe = Block22Probe(Block22Probe(:,2)>Block22ScanStart,:);

        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);

        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);

        Block22Probe(:,2) = Block22Probe(:,2)-Block22ScanStart;
        Block22Probe(:,3) = Block22Probe(:,3)-Block22ScanStart;
        Block22Probe(:,9) = Block22Probe(:,3)-Block22Probe(:,2);
        Block22Probe(:,8) = Block22Probe(:,8)-Block22ScanStart;
        Block22Probe(:,10) = Block22Probe(:,3)-Block22Probe(:,8);
        Block22Probe = array2table(Block22Probe,...
            'VariableNames',{'ProbeNum','ProbeOnset','ProbeEnd', 'PrevTrialNum', 'ProbeRes1', 'ProbeRes2', 'ProbeRes3','Q3Onset', 'ProbeDur', 'Q3Dur'});
    case (strcmp(C{3}, 's208')==1)
        Block2Diff=(secBlockStart-EEGTRIG.B1_list(2))/5000;
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(forthBlockStart-Block4_firstS1)/5000;

        Block2ScanStart = TaskLog.all_tstartblock(1,2)+Block2Diff;
        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;

        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)<3, [1,3,4,7,17:19, 16]);
        Block1Probe = Block1Probe(1:18,:);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block2Probe = Block2Probe(10,:);
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);

        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);
        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);

    case (strcmp (C{3}, 's214')==1)
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block2Diff=(secBlockStart-Block2_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(forthBlockStart-14258045)/5000;

        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block2ScanStart = TaskLog.all_tstartblock(1,2)-Block2Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;
        
        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==1, [1,3,4,7,17:19, 16]);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);
        Block4Probe = Block4Probe(2:end,:);

        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);
        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);

    case (strcmp(C{3}, 's218')==1)
        Block2Diff=(secBlockStart-EEGTRIG.B1_list(2))/5000;
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(forthBlockStart-Block4_firstS1)/5000;

        Block2ScanStart = TaskLog.all_tstartblock(1,2)+Block2Diff;
        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;

        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)<3, [1,3,4,7,17:19, 16]);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block1Probe = Block1Probe(Block1Probe(:,2)<scanEnd,:); 
        Block1Probe = Block1Probe(1:end-1,:); % remove unusable probe
        Block2Probe = Block2Probe(Block2Probe(:,2)>0,:); 
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);

        
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);
        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);

  case (strcmp(C{3}, 's220')==1)
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block2Diff=(secBlockStart-Block2_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(12966023-Block4_firstS1)/5000;
        
        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block2ScanStart = TaskLog.all_tstartblock(1,2)-Block2Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;
        
        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==1, [1,3,4,7,17:19, 16]);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);
    
        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);
        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);



   case (strcmp(C{3}, 's228')==1)

        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block2Diff=(secBlockStart-Block2_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(EEGTRIG.B1_list(4)-Block4_firstS1)/5000;
        
        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block2ScanStart = TaskLog.all_tstartblock(1,2)-Block2Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;

        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==1, [1,3,4,7,17:19, 16]);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)>=3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);

        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block3Probe = Block3Probe(Block3Probe(:,2)<scanEnd,:); 
        Block4Probe = Block4Probe(Block4Probe(:,2)>0,:); 
        
        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);

    otherwise
        Block1Diff=(firstBlockStart-Block1_firstS1)/5000;
        Block2Diff=(secBlockStart-Block2_firstS1)/5000;
        Block3Diff=(thirdBlockStart-Block3_firstS1)/5000;
        Block4Diff=(forthBlockStart-Block4_firstS1)/5000;
        
        Block1ScanStart = TaskLog.all_tstartblock(1,1)-Block1Diff;
        Block2ScanStart = TaskLog.all_tstartblock(1,2)-Block2Diff;
        Block3ScanStart = TaskLog.all_tstartblock(1,3)-Block3Diff;
        Block4ScanStart = TaskLog.all_tstartblock(1,4)-Block4Diff;
        
        Block1Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==1, [1,3,4,7,17:19, 16]);
        Block2Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==2, [1,3,4,7,17:19, 16]);
        Block3Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==3, [1,3,4,7,17:19, 16]);
        Block4Probe = TaskLog.probe_res(TaskLog.probe_res(:,5)==4, [1,3,4,7,17:19, 16]);
    
        Block1Probe(:,2) = Block1Probe(:,2)-Block1ScanStart;
        Block1Probe(:,3) = Block1Probe(:,3)-Block1ScanStart;
        Block1Probe(:,9) = Block1Probe(:,3)-Block1Probe(:,2);
        Block1Probe(:,8) = Block1Probe(:,8)-Block1ScanStart;
        Block1Probe(:,10) = Block1Probe(:,3)-Block1Probe(:,8);
        Block2Probe(:,2) = Block2Probe(:,2)-Block2ScanStart;
        Block2Probe(:,3) = Block2Probe(:,3)-Block2ScanStart;
        Block2Probe(:,8) = Block2Probe(:,8)-Block2ScanStart;
        Block2Probe(:,9) = Block2Probe(:,3)-Block2Probe(:,2);
        Block2Probe(:,10) = Block2Probe(:,3)-Block2Probe(:,8);
        Block3Probe(:,2) = Block3Probe(:,2)-Block3ScanStart;
        Block3Probe(:,3) = Block3Probe(:,3)-Block3ScanStart;
        Block3Probe(:,8) = Block3Probe(:,8)-Block3ScanStart;
        Block3Probe(:,9) = Block3Probe(:,3)-Block3Probe(:,2);
        Block3Probe(:,10) = Block3Probe(:,3)-Block3Probe(:,8);
        Block4Probe(:,2) = Block4Probe(:,2)-Block4ScanStart;
        Block4Probe(:,3) = Block4Probe(:,3)-Block4ScanStart;
        Block4Probe(:,8) = Block4Probe(:,8)-Block4ScanStart;
        Block4Probe(:,9) = Block4Probe(:,3)-Block4Probe(:,2);
        Block4Probe(:,10) = Block4Probe(:,3)-Block4Probe(:,8);
end

% Catch late scanning start
if any (Block1Probe(:,2)< 0)
   fprintf('Timing inconsitency in Probe Block 1 for %s \n',C{3});
   Block1Probe = Block1Probe(Block1Probe(:,2)>0, :);
end

if any (Block2Probe(:,2)< 0)
   fprintf('Timing inconsitency in Probe Block 2 for %s \n',C{3});
   Block2Probe = Block2Probe(Block2Probe(:,2)>0, :);
end

if any (Block3Probe(:,2)< 0)
   fprintf('Timing inconsitency in Probe Block 3 for %s \n',C{3});
   Block3Probe = Block3Probe(Block3Probe(:,2)>0, :);
end 

if any (Block4Probe(:,2)< 0)
   Block4Probe = Block4Probe(Block4Probe(:,2)>0, :);
end

Block1Probe = array2table(Block1Probe,...
    'VariableNames',{'ProbeNum','ProbeOnset','ProbeEnd', 'PrevTrialNum', 'ProbeRes1', 'ProbeRes2', 'ProbeRes3','Q3Onset', 'ProbeDur', 'Q3Dur'});
Block2Probe = array2table(Block2Probe,...
    'VariableNames',{'ProbeNum','ProbeOnset','ProbeEnd', 'PrevTrialNum', 'ProbeRes1', 'ProbeRes2', 'ProbeRes3', 'Q3Onset', 'ProbeDur', 'Q3Dur'});
Block3Probe = array2table(Block3Probe,...
    'VariableNames',{'ProbeNum','ProbeOnset','ProbeEnd', 'PrevTrialNum', 'ProbeRes1', 'ProbeRes2', 'ProbeRes3','Q3Onset', 'ProbeDur',  'Q3Dur'});
Block4Probe = array2table(Block4Probe,...
    'VariableNames',{'ProbeNum','ProbeOnset','ProbeEnd', 'PrevTrialNum', 'ProbeRes1', 'ProbeRes2', 'ProbeRes3','Q3Onset',  'ProbeDur', 'Q3Dur'});

%%

% %savepath = '/Volumes/TOSHIBA/MW_EEGFMRI_PROJECT/EEGFMRI_BehavData/TimingforfMRI/';
savepath =[path 'Shared drives/MW_fMRI_EEG/Data/MRI/TimingforfMRI/'];

% create folder
if exist([savepath C{3}]) == 0
    fprintf('Folder does not exist. Creating: %s \n',C{3})
    mkdir([savepath C{3}])
    mkdir([savepath C{3} '/Block1'])
    mkdir([savepath C{3} '/Block2'])
    mkdir([savepath C{3} '/Block3'])
    mkdir([savepath C{3} '/Block4'])
else
    fprintf('%s Folder exists. \n',C{3})
end

savefilenameB1 = ['behavTimes_', C{3} '_Block1.csv'];
savefilenameB2 = ['behavTimes_', C{3} '_Block2.csv'];
savefilenameB3 = ['behavTimes_', C{3} '_Block3.csv'];
savefilenameB4 = ['behavTimes_', C{3} '_Block4.csv'];
savefilenameB1P = ['behavTimes_', C{3} '_ProbeBlock1.csv'];
savefilenameB2P = ['behavTimes_', C{3} '_ProbeBlock2.csv'];
savefilenameB3P = ['behavTimes_', C{3} '_ProbeBlock3.csv'];
savefilenameB4P = ['behavTimes_', C{3} '_ProbeBlock4.csv'];


writetable(Block1,[[savepath C{3} '/Block1'] filesep savefilenameB1]);
writetable(Block2,[[savepath C{3} '/Block2'] filesep savefilenameB2]);
writetable(Block3,[[savepath C{3} '/Block3'] filesep savefilenameB3]);
writetable(Block4,[[savepath C{3} '/Block4'] filesep savefilenameB4]);
writetable(Block1Probe,[[savepath C{3} '/Block1'] filesep savefilenameB1P]);
writetable(Block2Probe,[[savepath C{3} '/Block2'] filesep savefilenameB2P]);
writetable(Block3Probe,[[savepath C{3} '/Block3'] filesep savefilenameB3P]);
writetable(Block4Probe,[[savepath C{3} '/Block4'] filesep savefilenameB4P]);

if  strcmp (C{3}, 's207')
    mkdir([savepath C{3} '/Block22'])
    savefilenameB22 = ['behavTimes_', C{3} '_Block22.csv'];
    savefilenameB22P = ['behavTimes_', C{3} '_ProbeBlock22.csv'];
    writetable(Block22,[[savepath C{3} '/Block22'] filesep savefilenameB22]);
    writetable(Block22Probe,[[savepath C{3} '/Block22'] filesep savefilenameB22P]);   
end

end

