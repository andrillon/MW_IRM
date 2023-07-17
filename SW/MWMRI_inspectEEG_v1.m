%%
run ../localdef.m
addpath((path_fieldtrip));
ft_defaults;

%% choose and load subject
List_Subj=dir([data_path filesep '*' filesep '*' filesep '*EO*.bdf']);
ListNames={List_Subj.name};
pick=listdlg('ListString',ListNames);
fprintf('... working on %s\n',ListNames{pick})

%% basic preprocessing

addpath(genpath(path_eeglab));
EEG = pop_loadset( 'filename',[files(nF).folder filesep files(nF).name]);
rmpath(genpath(path_eeglab));
eeglab2f    
    
hdr=ft_read_header([List_Subj(pick).folder filesep List_Subj(pick).name]);
cfg=[];
cfg.channel        = hdr.label(match_str(hdr.chantype,'eeg'));
cfg.dataset        = [List_Subj(pick).folder filesep List_Subj(pick).name];
cfg.demean         = 'yes';
cfg.lpfilter       = 'yes';        % enable high-pass filtering
cfg.lpfilttype     = 'but';
cfg.lpfiltord      = 5;
cfg.lpfreq         = 40;
cfg.hpfilter       = 'yes';         % enable high-pass filtering
cfg.hpfilttype     = 'but';
cfg.hpfiltord      = 5;
cfg.hpfreq         = 0.1;
cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
cfg.dftfreq        = [50]; % set up the frequencies for notch filtering
data               = ft_preprocessing(cfg); % read raw data

%% reject trials
cfg          = [];
cfg.method   = 'summary';
cfg.alim     = 5e-5;
data2        = ft_rejectvisual(cfg,data);

%% display data
cfg=[];
cfg.continuous='no';
cfg.allowoverlap='true';
cfg.viewmode='vertical';
cfg = ft_databrowser(cfg, data);
