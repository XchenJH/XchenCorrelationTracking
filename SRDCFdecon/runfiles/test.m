function results=test(seq, res_path, bSaveImage, parameters)

close all

s_frames = seq.s_frames;

%translation filter feature parameters
hog_params.nDim = 31;

grayscale_params.colorspace='gray';
grayscale_params.nDim = 1;

cn_params.tablename = 'CNnorm';
cn_params.nDim = 10;
cn_params.useForGray = false;

ic_params.tablename = 'intensityChannelNorm6';
ic_params.nDim = 5;
ic_params.useForColor = false;

params.t_features = {
    ...struct('getFeature',@get_colorspace, 'fparams',grayscale_params),...
    struct('getFeature',@get_fhog,'fparams',hog_params),...
    ...struct('getFeature',@get_table_feature, 'fparams',cn_params),...
    ...struct('getFeature',@get_table_feature, 'fparams',ic_params),...
};
params.t_global.cell_size = 4; %global feature-cell size
params.t_global.cell_selection_thresh = 0.75^2;

% Translation filter parameters
params.search_area_shape = 'square';    % Shape of the training region
params.search_area_scale = 4.0;         % the scaling of the target size to get the search area
params.output_sigma_factor = 1/16;		% spatial bandwidth (proportional to target)
params.refinement_iterations = 1;       % number of iterations used to refine the resulting position in a frame
params.translation_model_max_area = 50^2;   % the maximal search area in number of feature grid cells
params.interpolate_response = 4;        % correlation score interpolation strategy
params.gradient_ascent_iterations = 5;  % Number of Newton iteration
params.init_strategy = 'indep';         % strategy for initializing the filter, 'const_reg' or 'indep'
params.num_init_GS_iter = 4;            % Number of initial Gauss-Seidel iterations
params.num_GS_iter = 4;                 % Number of Gauss-Seidel iterations
params.num_acs_iter = 1;                % Number of Alternate Convex Search iterations

% Regularization window parameters
params.use_reg_window = 1;              % wather to use windowed regularization or not
params.reg_window_min = 0.1;			% the minimum value of the regularization window
params.reg_window_edge = 3.0;           % the impact of the spatial regularization (value at the target border), depends on the detection size and the feature dimensionality
params.reg_window_power = 2;            % the degree of the polynomial to use (e.g. 2 is a quadratic window)
params.reg_sparsity_threshold = 0.05;   % a relative threshold of which DFT coefficients that should be set to zero
params.lambda = 1e-2;					% the weight of the standard (uniform) regularization, only used when params.use_reg_window == 0

% Scale parameters
params.number_of_scales = 7;
params.scale_step = 1.01;

% Learning parameters
params.learning_rate = 0.035;			% learning rate
params.nSamples = 300;                  % Maximal number of samples in memory
params.sample_reg = 5.0;                % Weights regularization (mu)
params.sample_burnin = 10;              % Number of frames before weight optimization starts
params.sample_weight_after_filter = true; % Do the weight optimization before or after filter optimization
params.sample_replace_strategy = 'constant_tail';   % How to set the prior weights
params.lt_size = 250;                   % Number of samples in the "constant_tail"

% Other parameters
params.visualization = 1;
params.debug = 1;

params.visualization = params.visualization || params.debug;

params.wsize = [seq.init_rect(1,4), seq.init_rect(1,3)];
params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);
params.s_frames = s_frames;

results = tracker(params);
