% Run this file to build the needed mex-files on windows

[home_dir, name, ext] = fileparts(mfilename('fullpath'));
cd(home_dir)

% Build merResize
cd mexResize
mex -lopencv_core242 -lopencv_imgproc242 -L./ -I./ mexResize.cpp MxArray.cpp
cd(home_dir)

% Build setnonzeros from the lightspeed matlab toolbox
cd external_libs
flags = ' -largeArrayDims ';
% flags = ' -DmwSize=int -DmwIndex=int '; % Use for really old versions of Matlab
eval(['mex' flags 'setnonzeros.c'])
cd(home_dir)

% Build gradientMex from Piotrs toolbox
cd feature_extraction
mex gradientMex.cpp -I./
cd(home_dir)

% Build mtimesx
cd external_libs
mtimesx_build
cd(home_dir)