

%% Init
% clear all;
% close all;

run /Users/kuszti/Documents/Git/MW_IRM/SW/localdef.m
files=dir([save_path filesep 'prct_RPA_SW_*']);

for nF=21 %1:length(files)

saveWaves =load([files(nF).folder filesep files(nF).name]);
slow_Waves= saveWaves.slow_Waves;


%%
% Change onset and negative peak time into miliseconds from the probe
slow_Waves(:,5)=slow_Waves(:,5)-(25*500);
slow_Waves(:,8)=slow_Waves(:,8)-(25*500);
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
BlockDataset = dir([TimingPath subFolder '/*/*_Block*.csv']); 
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
Block1sws(:,16)=Block1sws(:,2);
Block2sws(:,16)=Block2sws(:,2)-10;
Block3sws(:,16)=Block3sws(:,2)-20;
Block4sws(:,16)=Block4sws(:,2)-30;

% this function 
[Block1SWF, Block1SWC, Block1SWB] = SW_TimeRecodingforMR (Block1sws, BlockP{1},saveWaves);
[Block2SWF, Block2SWC, Block2SWB] = SW_TimeRecodingforMR (Block2sws, BlockP{2},saveWaves);
[Block3SWF, Block3SWC, Block3SWB] = SW_TimeRecodingforMR (Block3sws, BlockP{3},saveWaves);

if strcmp(subFolder, 's228')==1
    fprintf('No SWs in Block4 in subject: %s \n',subFolder);
else
    [Block4SWF, Block4SWC, Block4SWB] = SW_TimeRecodingforMR (Block4sws, BlockP{4},saveWaves);
end
%%
% SAVE

% BLOCK 1
savefilenameB1SF = ['swTimes_', subFolder '_SWBlock1F.csv'];
writetable(Block1SWF,[[TimingPath subFolder '/Block1'] filesep savefilenameB1SF]);
savefilenameB1SC = ['swTimes_', subFolder '_SWBlock1B.csv'];
writetable(Block1SWC,[[TimingPath subFolder '/Block1'] filesep savefilenameB1SC]);
savefilenameB1SB = ['swTimes_', subFolder '_SWBlock1C.csv'];
writetable(Block1SWB,[[TimingPath subFolder '/Block1'] filesep savefilenameB1SB]);

% BLOCK 2
savefilenameB2SF = ['swTimes_', subFolder '_SWBlock2F.csv'];
writetable(Block2SWF,[[TimingPath subFolder '/Block2'] filesep savefilenameB2SF]);
savefilenameB2SC = ['swTimes_', subFolder '_SWBlock2B.csv'];
writetable(Block2SWC,[[TimingPath subFolder '/Block2'] filesep savefilenameB2SC]);
savefilenameB2SB = ['swTimes_', subFolder '_SWBlock2C.csv'];
writetable(Block2SWB,[[TimingPath subFolder '/Block2'] filesep savefilenameB2SB]);


% BLOCK 3
savefilenameB3SF = ['swTimes_', subFolder '_SWBlock3F.csv'];
writetable(Block3SWF,[[TimingPath subFolder '/Block3'] filesep savefilenameB3SF]);
savefilenameB3SC = ['swTimes_', subFolder '_SWBlock3B.csv'];
writetable(Block3SWC,[[TimingPath subFolder '/Block3'] filesep savefilenameB3SC]);
savefilenameB3SB = ['swTimes_', subFolder '_SWBlock3C.csv'];
writetable(Block3SWB,[[TimingPath subFolder '/Block3'] filesep savefilenameB3SB]);




if strcmp(subFolder, 's228')==1

else
    % BLOCK 4
savefilenameB4SF = ['swTimes_', subFolder '_SWBlock4F.csv'];
writetable(Block4SWF,[[TimingPath subFolder '/Block4'] filesep savefilenameB4SF]);
savefilenameB4SC = ['swTimes_', subFolder '_SWBlock4B.csv'];
writetable(Block4SWC,[[TimingPath subFolder '/Block4'] filesep savefilenameB4SC]);
savefilenameB4SB = ['swTimes_', subFolder '_SWBlock4C.csv'];
writetable(Block4SWB,[[TimingPath subFolder '/Block4'] filesep savefilenameB4SB]);
end

end





