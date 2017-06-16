
% This demo script runs the SRDCFdecon tracker on the included "Couple" video.

% Add paths
setup_paths();

% Load video information
video_path = 'sequences/Skating2';
% video_path = 'sequences/Couple';
[seq, ~] = load_video_info(video_path);

% Run the tracker using the test.m runfile
results = test(seq);