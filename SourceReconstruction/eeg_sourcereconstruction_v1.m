%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Source Reconstruction for EEG/fMRI Mind-wandering project
%
% Author: Aniko Kusztor (aniko.kusztor@monash.edu)
%
% Last Update: 2024-01-25
%
% Created: 2024-01-2, Melbourne, Australia
% Version: MATLAB '23.2.0.2485118 (R2023b) Update 6'
% Operating System: macOS  Version: 14.2.1 Build: 23C71 
% Java Version: Java 1.8.0_402-b08 with Amazon.com Inc. OpenJDK 64-Bit Server VM mixed mode
% FieldTrip version is: 20240110
% FieldTrip-SimBio pipeline: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% mex -setup -v


run ../localdef.m

% adding relevant toolboxes to the path
addpath((path_fieldtrip))
ft_defaults
%addpath(genpath(path_eeglab));
addpath((ft_spm12_path));
addpath(genpath(freesurfer_path));
addpath(genpath(simbio_path));


MRIpath='/Volumes/MNHS-SPP/EEG-fMRI-MW/mindwander_eeg/Data/Analysis/s217/anat/sMRH119_MWMRI217_MR01-0002-00001-000192-01.nii';


% HEADMODEL
% read mri
mri = ft_read_mri(MRIpath);

cfg = [];
cfg.method = 'interactive';
mri = ft_volumerealign(cfg, mri);


% Using the n/l/r keys the fiducials can be specified, the
% landmarks can be specified with a/p/z. When pressing q the interactive
% mode will stop and the transformation matrix is computed.
% n -- Nasion --
% l & r --Left/Right Pre Auricular point (LPA/RPA)--
% a --
% p --
% z -- top of the skull

% reslice
cfg = [];
cfg.method = 'flip';
mri = ft_volumereslice(cfg,mri);

%save mri_resliced.mat mri

% segmentation
cfg = [];
cfg.output = {'gray','white','csf','skull','scalp'};
cfg.spmversion = 'spm12';
cfg.scalpthreshold = 0.06;
segmentedmri = ft_volumesegment(cfg, mri);

%save segmentedmri segmentedmri

% plot segmentation image
seg_i            = ft_datatype_segmentation(segmentedmri,'segmentationstyle',...
                                            'indexed');
cfg              = [];
cfg.funparameter = 'seg';
cfg.location     = 'center';

% % convert from probabilistic/binary into indexed representation
% segmentedmri_indexed = ft_datatype_segmentation(segmentedmri, 'segmentationstyle', 'indexed');
% 
% % also add the anatomical mri
% segmentedmri_indexed.anatomy = mri_resliced.anatomy;
% 
% cfg              = [];
% cfg.anaparameter = 'anatomy';
% cfg.funparameter = 'tissue';
% cfg.funcolormap  = lines(6);              % distinct color per tissue + background
% cfg.atlas        = segmentedmri_indexed;  % this is just like an anatomical atlas, see https://www.fieldtriptoolbox.org/template/atlas/
% cfg.location     = 'center';
% ft_sourceplot(cfg, segmentedmri_indexed);


% create mesh
cfg=[];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg,segmentedmri);

save mesh.mat mesh

disp(mesh)

% create headmodel
cfg = [];
cfg.method = 'simbio';
cfg.conductivity = [0.33 0.14 1.79 0.01 0.43]; % same as tissuelabel in vol
vol = ft_prepare_headmodel(cfg, mesh);

save vol.mat vol

disp(vol)


%%
% Brainstorm pipeline
cd /Applications/brainstorm3/
brainstorm

addpath('/Users/kuszti/Applications/SimNIBS-4.0/bin')


%setenv('LD_LIBRARY_PATH', sprintf('/Users/kuszti/Applications/SimNIBS-4.0/simnibs_env/lib/python3.9/site-packages/simnibs/external/lib/osx/', getenv('LD_LIBRARY_PATH')))


setenv('PATH', [getenv('PATH') '/Users/kuszti/Applications/SimNIBS-4.0/bin']);


