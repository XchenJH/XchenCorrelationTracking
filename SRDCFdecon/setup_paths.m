function setup_paths()

% Add the neccesary paths

[pathstr, name, ext] = fileparts(mfilename('fullpath'));

% Tracker implementation
addpath([pathstr '/implementation/']);

% Runfiles
addpath([pathstr '/runfiles/']);

% Utilities
addpath([pathstr '/utils/']);

% The feature extraction
addpath(genpath([pathstr '/feature_extraction/']));

% Mtimesx
addpath([pathstr '/external_libs/']);

% mexResize
addpath([pathstr '/mexResize/']);