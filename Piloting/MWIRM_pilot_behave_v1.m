%%
clear all
close all

run ../localdef
addpath(genpath(lscpTools_path))

% look for behav files
data_path=[root_path filesep 'Pilot/'];
files=dir([data_path filesep 'wanderIM_behavres_s1*.mat']);

%%
all_probe_res=[];
all_test_res=[];
for n=1:length(files)
    % load
    load([data_path filesep files(n).name]);
    SubID=SubjectInfo.subID;
    fprintf('... %s duration of task %3.1f min\n',SubID,(all_GrandEnd-all_GrandStart)/60)
    
    % Performance on SART
    %  1: num block
    %  2: block cond 1: digit / 2: emotions / 3: gender
    %  3: image set
    %  4: num trial
    %  5: seq trial
    %  6: target
    %  7: resp
    %  8: stim onset
    %  9: stim pre
    % 10: resp onset
    % 11: nogo
    % 12: go
    RTs=[test_res(:,10)-test_res(:,8)];
    
    tp_nogos=test_res(~isnan(test_res(:,11)),11);
    tp_gos=test_res(~isnan(test_res(:,12)),12);
    [dprime_test, crit_test]=calc_dprime((tp_gos==1),(tp_nogos==0));
    corr_go=nanmean(tp_gos);
    corr_nogo=nanmean(tp_nogos);
    
    rt_gos=nanmean(RTs(~isnan(test_res(:,12))));
    rt_nogos=nanmean(RTs(~isnan(test_res(:,11))));
    all_test_res=[all_test_res ; [n unique(test_res(:,2)) corr_go corr_nogo dprime_test crit_test rt_gos rt_nogos]];
    
    % Data on probe
    %  1: num probe within block
    %  2: probe onset (expected)
    %  3: probe onset (actual)
    %  4: num block
    %  5: block type (1: digit / 2: emotions / 3: gender)
    %  6: num trial
    %  7: resp Q1 Mind-State (1: On / 2: MW / 3: MB / 4: ?)
    %  8: resp Q2 Origin (1: room / 2: personal / 3: task)
    %  9: resp Q3 Alterness (1 (alert) to 4 (sleepy))
    % 10-12: onset R1-3
    % 13-15: onset Q1-3
    % 16-18: resp number Q1-3
    temp=probe_res(:,[5 16:18]);
    clear pption pption2 vigilance
    for nstate=1:4
        pption(nstate)=mean(temp(:,2)==nstate);
        if nstate~=4
            pption2(nstate)=mean(temp(temp(:,3)==2,3)==nstate);
        end
        vigilance(nstate)=nanmean(temp(temp(:,2)==nstate,4));
    end
    pption_MS=pption;
    pption_Origin=pption2;
    vigilance_MS=vigilance;
    all_probe_res=[all_probe_res ; [n unique(probe_res(:,5)) pption_MS pption_Origin vigilance_MS]];
end

%% COnvert to table
%  all_test_res=[all_test_res ; [n unique(test_res(:,2)) corr_go corr_nogo dprime_test crit_test rt_gos rt_nogos]];
%  all_probe_res=[all_probe_res ; [n unique(probe_res(:,5)) pption_MS pption_MS vigilance_MS]];
test_table=array2table(all_test_res,'VariableNames',{'nFile','Task','GoPerf','NoGoPerf','dp','crit','RTgo','RTnogo'});
test_table.nFile=categorical(test_table.nFile);
test_table.Task=categorical(test_table.Task);
test_table.Task(test_table.Task=='1')='Digit';
test_table.Task(test_table.Task=='2')='Emotion';
test_table.Task(test_table.Task=='3')='Gender';


probe_table=array2table(all_probe_res,'VariableNames',{'nFile','Task','ppON','ppOFF','ppMB','ppDR','ppDIS','ppMW','ppINT','vigON','vigOFF','vigMB','vigDR'});
probe_table.nFile=categorical(probe_table.nFile);
probe_table.Task=categorical(probe_table.Task);
probe_table.Task(probe_table.Task=='1')='Digit';
probe_table.Task(probe_table.Task=='2')='Emotion';
probe_table.Task(probe_table.Task=='3')='Gender';

%% Plots
TaskType={'Digit','Emotion','Gender'};
figure;
for nplot=1:6
    subplot(3,2,nplot); hold on;
    switch nplot
        case 1
            temp=test_table.crit;
            title('criterion')
        case 2
            temp=test_table.dp;
            title('d''')
        case 3
            temp=test_table.GoPerf;
            title('Go Perf')
        case 4
            temp=test_table.NoGoPerf;
            title('NoGo Perf')
        case 5
            temp=test_table.RTgo;
            title('RT Go')
        case 6
            temp=test_table.RTnogo;
            title('RT NoGo')
    end
    temp2=[];
    for ntask=1:3
        scatter(ntask,temp(test_table.Task==TaskType{ntask}),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','SizeData',144)
        temp2(ntask)=mean(temp(test_table.Task==TaskType{ntask}));
    end
    plot(1:3,temp2,'LineWidth',3,'Color','k');
    set(gca,'XTick',1:3,'XTickLabel',TaskType); xlim([0.5 3.5]); format_fig;
end

figure;
for nplot=1:2
    subplot(1,2,nplot); hold on;
    switch nplot
        case 1
            startidx=2;
            title('PPtion Mind-State')
        case 2
            startidx=9;
            title('Tiredness')
    end
    hp=[];
    for nstate=1:3
        temp=table2array(probe_table(:,startidx+nstate));
        temp2=[];
        for ntask=1:3
            scatter(ntask,temp(probe_table.Task==TaskType{ntask}),'Marker','o','MarkerEdgeColor',Colors(nstate,:),'MarkerFaceColor',Colors(nstate,:),'SizeData',144)
            temp2(ntask)=mean(temp(probe_table.Task==TaskType{ntask}));
        end
        hp(nstate)=plot(1:3,temp2,'LineWidth',3,'Color',Colors(nstate,:));
    end
    set(gca,'XTick',1:3,'XTickLabel',TaskType); xlim([0.5 3.5]); format_fig;
end
legend(hp,{'ON','OFF','MB'});