path_ROI='/Users/thandrillon/Data/MW_IRM/MRI/ROI_timecourse/';
path_Behav='/Users/thandrillon/Data/MW_IRM/Behav/';
files=dir([path_ROI filesep 'rffa_epoched30s10s_*.mat']);

all_ProbeTrials=[];
table_ROI=array2table(zeros(0,9),'VariableNames',{'SubID','Block','Probe','totProbe','MS','Vig','ROI_w30','ROI_w20','ROI_w10'});
table_ROI.SubID=categorical(table_ROI.SubID);
table_ROI.MS=categorical(table_ROI.MS);

xTime=-30:1.5769:10;

%%% Summary of subject information: Behavioural/fMRI/EEG Data Summary
%%% https://docs.google.com/spreadsheets/d/14yzZas4tCX-DoXkF9rvfsoEVWjokEmuqUHZkgWFtZPo/edit#gid=0

for nF=1:length(files)
    load([path_ROI filesep files(nF).name])
    split_name=split(files(nF).name,'_');
    
    behav_file=dir([path_Behav filesep 'wanderIM_behavres_' split_name{end}(1:end-4) '*.mat']);
    if ~isempty(behav_file)
        load([path_Behav filesep behav_file(end).name])
        if size(probe_res,1)==size(ProbeTrials,1)
            for nP=1:size(probe_res,1)
                size_table=size(table_ROI,1);
                table_ROI.SubID(size_table+1)=split_name{end}(1:end-4);
                table_ROI.Block(size_table+1)=probe_res(nP,5);
                table_ROI.Probe(size_table+1)=probe_res(nP,1);
                table_ROI.totProbe(size_table+1)=nP;
                if probe_res(nP,17)==1
                    table_ROI.MS(size_table+1)='ON';
                elseif probe_res(nP,17)==2
                    table_ROI.MS(size_table+1)='MW';
                elseif probe_res(nP,17)==3 || probe_res(nP,17)==4
                    table_ROI.MS(size_table+1)='MB';
                end
                
                table_ROI.Vig(size_table+1)=probe_res(nP,19);
                table_ROI.ROI_win30(size_table+1)=nanmean(ProbeTrials(nP,xTime>-30 & xTime<-20));
                table_ROI.ROI_win20(size_table+1)=nanmean(ProbeTrials(nP,xTime>-20 & xTime<-10));
                table_ROI.ROI_win10(size_table+1)=nanmean(ProbeTrials(nP,xTime>-10 & xTime<0));
            end
            %             all_ProbeTrials=cat(1,all_ProbeTrials,ProbeTrials);
        else
                    warning('Wrong number of lines')
        end
    else
        warning('Could not find behav file')
    end
    all_ProbeTrials=cat(1,all_ProbeTrials,ProbeTrials);
end