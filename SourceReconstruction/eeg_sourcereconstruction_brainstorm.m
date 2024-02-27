%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Source Reconstruction for EEG/fMRI Mind-wandering project
%
% Author: Aniko Kusztor (aniko.kusztor@monash.edu)
%
% Last Update: 2024-01-25
%
% Created: 2024-01-2, Melbourne, Australia
% Version: MATLAB 9.14.0.2489007 (R2023a) Update 6
% Operating System: macOS  Version: 14.2.1 Build: 23C71 
% Java Version: Java 1.8.0_402-b08 with Amazon.com Inc. OpenJDK 64-Bit Server VM mixed mode
% Brainstorm version is: 30-Jan-2024
% Brainstorm:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
close all
data_path='/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/EEG/Preprocessed with BVA';
MRIpath='/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/MRI/anat/';
reports_dir = '/Users/kuszti/Library/CloudStorage/GoogleDrive-aniko.kusztor@monash.edu/Shared drives/MW_fMRI_EEG/Data/';

% create file list as cells
EEGfileslist = dir([data_path filesep filesep 'MWMRI*clean5.set']);
MRIfileslist = dir([MRIpath filesep '*/anat/s*.nii']);
EEGfiles= fullfile({EEGfileslist.folder}, {EEGfileslist.name})';
MRIfiles= fullfile({MRIfileslist.folder}, {MRIfileslist.name})';


% Brainstorm pipeline
cd /Applications/brainstorm3/
%brainstorm
% ===== CREATE PROTOCOL =====
% The protocol name has to be a valid folder name (no spaces, no weird characters...)
ProtocolName = 'MWMRI_full';
% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end

% select subjects based on folder names
Subjects = dir(MRIpath);
Subjects = Subjects(4:end,:); %remove '.' '..' 'DS_store'

% % Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);


for nF=1:length(Subjects)
fprintf('... subject: %s\n',Subjects(nF).name)

% Input files
sFiles = [];
% Start a new report
bst_report('Start', sFiles);
% ===== IMPORT ANATOMY =====
% ===== IMPORT MRI VOLUMES =====
% Create subject
[sSubject, iSubject] = db_add_subject(Subjects(nF).name, [], 0, 0);
% Import T1 MRI
% Process: Import MRI
sFiles = bst_process('CallProcess', 'process_import_mri', sFiles, [], ...
    'subjectname', Subjects(nF).name, ...
    'mrifile',     {MRIfiles{nF}, 'Nifti1'}, ...
    'nas',         [0, 0, 0], ...
    'lpa',         [0, 0, 0], ...
    'rpa',         [0, 0, 0], ...
    'ac',          [0, 0, 0], ...
    'pc',          [0, 0, 0], ...
    'ih',          [0, 0, 0]);
% set anatomical landmarks ? 

% Process: Segment MRI with CAT12
sFiles = bst_process('CallProcess', 'process_segment_cat12', sFiles, [], ...
    'subjectname', Subjects(nF).name, ...
    'nvertices',   15000, ...
    'tpmnii',      {'', 'Nifti1'}, ...
    'sphreg',      1, ...
    'vol',         1, ...
    'extramaps',   1, ...
    'cerebellum',  0);

% Process: MNI normalization
sFiles = bst_process('CallProcess', 'process_mni_normalize', sFiles, [], ...
    'subjectname', Subjects(nF).name, ...
    'method',      'maff8', ...  % maff8:Affine registration using SPM mutual information algorithm.Estimates a simple 4x4 linear transformation to the MNI space.Included in Brainstorm.
    'uset2',       0);

% import EEG data
% Process: Create link to raw file
sFilesRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
    'subjectname',    Subjects(nF).name, ...
    'datafile',       {EEGfiles{nF}, 'EEG-EEGLAB'}, ...
    'channelreplace', 1, ...
    'channelalign',   1, ...
    'evtmode',        'value');

% Process: Refine registration
sFilesRaw = bst_process('CallProcess', 'process_headpoints_refine', sFilesRaw, [], ...
    'tolerance', 0);

% Process: Project electrodes on scalp
bst_process('CallProcess', 'process_channel_project', sFilesRaw, []);

% Process: Snapshot: Sensors/MRI registration
bst_process('CallProcess', 'process_snapshot', sFilesRaw, [], ...
    'target',   1, ...  % Sensors/MRI registration
    'modality', 4, ...  % EEG
    'orient',   1, ...  % left
    'Comment',  'MEG/MRI Registration');

% FACE trials
% T (-200ms - 600ms)
% Process: Import MEG/EEG: Events
sFilesEpoch = bst_process('CallProcess', 'process_import_data_event', sFilesRaw, [], ...
    'subjectname',   Subjects(nF).name, ...
    'condition',     'Face', ...
    'eventname',     'T  1', ...
    'timewindow',    [], ...
    'epochtime',     [-0.2, 0.6], ...
    'split',         0, ...
    'createcond',    1, ...
    'ignoreshort',   1, ...
    'channelalign',  1, ...
    'usectfcomp',    1, ...
    'usessp',        1, ...
    'freq',          [], ...
    'baseline',      [], ...
    'blsensortypes', 'MEG, EEG');

% Process: Refine registration
sFilesEpoch = bst_process('CallProcess', 'process_headpoints_refine', sFilesEpoch, [], ...
    'tolerance', 0);

% Process: Project electrodes on scalp
bst_process('CallProcess', 'process_channel_project', sFilesEpoch, []);

% Process: Snapshot: Sensors/MRI registration
bst_process('CallProcess', 'process_snapshot', sFilesEpoch, [], ...
    'target',   1, ...  % Sensors/MRI registration
    'modality', 4, ...  % EEG
    'orient',   1, ...  % left
    'Comment',  'MEG/MRI Registration');

% Process: Average: By folder (subject average)
sFilesAVR = bst_process('CallProcess', 'process_average', sFilesEpoch, [], ...
    'avgtype',       3, ...  % By folder (subject average)
    'avg_func',      1, ...  % Arithmetic average:  mean(x)
    'weighted',      0, ...
    'keepevents',    0);


% for SW 
% add EEG markers based on 





% ===== SOURCE ANALYSIS: SURFACE =====

% Process: Generate BEM surfaces
bst_process('CallProcess', 'process_generate_bem', [], [], ...
    'subjectname', Subjects(nF).name, ...
    'nscalp',      1922, ...
    'nouter',      1922, ...
    'ninner',      1922, ...
    'thickness',   4, ...
    'method',      'brainstorm');  % Brainstorm


% Process: Compute head model
bst_process('CallProcess', 'process_headmodel', sFilesAVR, [], ...
    'sourcespace', 1, ...  % Cortex surface
    'eeg',         3, ...  % OpenMEEG BEM
    'openmeeg',    struct(...
         'BemSelect',    [0, 0, 1], ...
         'BemCond',      [1, 0.0125, 1], ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemFiles',     {{}}, ...
         'isAdjoint',    1, ...
         'isAdaptative', 1, ...
         'isSplit',      0, ...
         'SplitLength',  4000));

% Process: Compute noise covariance FOR FACE baseline
% Process: Compute covariance (noise or data)
sFiles = bst_process('CallProcess', 'process_noisecov', sFilesEpoch, [], ...
    'baseline',       [-0.2, -0.002], ...
    'datatimewindow', [-7.21644966e-16, 0.6], ...
    'sensortypes',    '', ...
    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       0, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    1);  % Replace


% bst_process('CallProcess', 'process_noisecov', sFilesRaw, [], ...
%     'baseline',       [-25, -0.002], ...
%     'datatimewindow', [-25, -0.002], ...
%     'sensortypes',    'EEG', ...
%     'target',         1, ...  % Noise covariance     (covariance over baseline time window)
%     'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
%     'identity',       0, ...
%     'copycond',       1, ...
%     'copysubj',       0, ...
%     'replacefile',    1);  % Replace

% Process: Compute sources [2018]
sAvgSrc = bst_process('CallProcess', 'process_inverse_2018', sFilesAVR, [], ...
    'output',  1, ...  % Kernel only: shared
    'inverse', struct(...
         'Comment',        'sLORETA: EEG', ...
         'InverseMethod',  'minnorm', ...
         'InverseMeasure', 'sloreta', ...
         'SourceOrient',   {{'free'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       0, ...
         'WeightExp',      0.5, ...
         'WeightLimit',    10, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       0.1, ...
         'SnrMethod',      'fixed', ...
         'SnrRms',         1e-06, ...
         'SnrFixed',       3, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'EEG'}}));

% Save and display report
% Save and display report
ReportFile = bst_report('Save', []);
bst_report('Export', ReportFile, reports_dir);
bst_report('Open', ReportFile);

end


%%



% Process: Snapshot: Sources (one time)
bst_process('CallProcess', 'process_snapshot', sAvgSrc, [], ...
    'target',    8, ...  % Sources (one time)
    'modality',  1, ...  % MEG (All)
    'orient',    3, ...  % top
    'time',      0, ...
    'threshold', 60, ...
    'Comment',   'Average Face');


% ===== SOURCE ANALYSIS: VOLUME =====

% Process: Compute head model
bst_process('CallProcess', 'process_headmodel', sFilesAVR, [], ...
    'Comment',     '', ...
    'sourcespace', 2, ...  % MRI volume
    'volumegrid',  struct(...
         'Method',        'isotropic', ...
         'nLayers',       17, ...
         'Reduction',     3, ...
         'nVerticesInit', 4000, ...
         'Resolution',    0.005, ...
         'FileName',      []), ...
    'eeg',         3, ...  % OpenMEEG BEM
    'openmeeg',    struct(...
         'BemSelect',    [1, 1, 1], ...
         'BemCond',      [1, 0.0125, 1], ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemFiles',     {{}}, ...
         'isAdjoint',    0, ...
         'isAdaptative', 1, ...
         'isSplit',      0, ...
         'SplitLength',  4000), ...
    'channelfile', '');

% Process: Compute sources [2018]
sAvgSrcVol = bst_process('CallProcess', 'process_inverse_2018', sFilesAVR, [], ...
    'output',  1, ...  % Kernel only: shared
    'inverse', struct(...
         'Comment',        'Dipoles: EEG', ...
         'InverseMethod',  'gls', ...
         'InverseMeasure', 'performance', ...
         'SourceOrient',   {{'free'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      0.5, ...
         'WeightLimit',    10, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       0.1, ...
         'SnrMethod',      'rms', ...
         'SnrRms',         1e-06, ...
         'SnrFixed',       3, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'EEG'}}));


%%
% Process: Snapshot: Sources (one time)
bst_process('CallProcess', 'process_snapshot', sAvgSrcVol, [], ...
    'target',    8, ...  % Sources (one time)
    'modality',  1, ...  % MEG (All)
    'orient',    3, ...  % top
    'time',      0, ...
    'threshold', 0, ...
    'Comment',   'Dipole modeling');

% Process: Dipole scanning
sDipScan = bst_process('CallProcess', 'process_dipole_scanning', sAvgSrcVol, [], ...
    'timewindow', [-0.040, 0.100], ...
    'scouts',     {});
% Process: Snapshot: Dipoles
bst_process('CallProcess', 'process_snapshot', sDipScan, [], ...
    'target',    13, ...  % Dipoles
    'orient',    3, ...   % top
    'threshold', 90, ...
    'Comment',   'Dipole scanning');
% Process: Snapshot: Dipoles
bst_process('CallProcess', 'process_snapshot', sDipScan, [], ...
    'target',    13, ...  % Dipoles
    'orient',    1, ...   % left
    'threshold', 90, ...
    'Comment',   'Dipole scanning');

% Process: FieldTrip: ft_dipolefitting
sDipFit = bst_process('CallProcess', 'process_ft_dipolefitting', sFilesAvg, [], ...
    'timewindow',  [-0.040, 0.100], ...
    'sensortypes', 'EEG', ...
    'dipolemodel', 1, ...  % Moving dipole
    'numdipoles',  1, ...
    'volumegrid',  [], ...
    'symmetry',    0, ...
    'filetag',     '');
% Process: Snapshot: Dipoles
bst_process('CallProcess', 'process_snapshot', sDipFit, [], ...
    'target',    13, ...  % Dipoles
    'orient',    3, ...   % top
    'threshold', 95, ...
    'Comment',   'Dipole fitting');
% Process: Snapshot: Dipoles
bst_process('CallProcess', 'process_snapshot', sDipFit, [], ...
    'target',    13, ...  % Dipoles
    'orient',    1, ...   % left
    'threshold', 95, ...
    'Comment',   'Dipole fitting');

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);


