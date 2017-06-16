function results = tracker(params)

search_area_scale = params.search_area_scale;                         	%extra area surrounding the target
output_sigma_factor = params.output_sigma_factor;	%spatial bandwidth (proportional to target)
lambda = params.lambda;                 			%regularization
learning_rate = params.learning_rate;       		%linear interpolation factor for adaptation
refinement_iterations = params.refinement_iterations;       %number of iterations used to refine the resulting position in a frame
translation_model_max_area = params.translation_model_max_area;
nScales = params.number_of_scales;
scale_step = params.scale_step;
interpolate_response = params.interpolate_response;
params.visualization = params.visualization || params.debug;

if ~isfield(params, 'num_acs_iter')
    params.num_acs_iter = 1;
end

features = params.t_features;

s_frames = params.s_frames;
pos = floor(params.init_pos);
target_sz = floor(params.wsize);

debug = params.debug;

num_frames = numel(s_frames);

init_target_sz = target_sz;

%set the feature ratio to the feature-cell size
featureRatio = params.t_global.cell_size;

search_area = prod(init_target_sz / featureRatio * search_area_scale);

if isfield(params.t_global, 'cell_selection_thresh')
    if search_area < params.t_global.cell_selection_thresh * translation_model_max_area
        params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * translation_model_max_area)))));
        
        featureRatio = params.t_global.cell_size;
        search_area = prod(init_target_sz / featureRatio * search_area_scale);
    end
end

global_feat_params = params.t_global;

if search_area > translation_model_max_area
    currentScaleFactor = sqrt(search_area / translation_model_max_area);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

%window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    case 'custom'
        sz = [base_target_sz(1)*2 base_target_sz(2)*4]; % for testing
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);

% construct the label function
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs] = ndgrid( rg,cg);
y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf = fft2(y);

if interpolate_response == 1
    interp_sz = use_sz * featureRatio;
else
    interp_sz = use_sz;
end

% construct cosine window
cos_window = single(hann(use_sz(1))*hann(use_sz(2))');

% the search area size
support_sz = prod(use_sz);

% Calculate feature dimension
im = imread(s_frames{1});
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end

feature_dim = 0;
for n = 1:length(features)
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end;
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end;
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        feature_dim = feature_dim + features{n}.fparams.nDim;
    end;
end;

if size(im,3) > 1 && colorImage == false
    im = im(:,:,1);
end

% compute the indices for the real, positive and negative parts of the
% spectrum
[dft_sym_ind, dft_pos_ind, dft_neg_ind] = partition_spectrum2(use_sz);
dft_nonneg_ind = [dft_sym_ind; dft_pos_ind];

% the discrete fourier series output indices
dfs_sym_ind = (1:length(dft_sym_ind))';
dfs_real_ind = dfs_sym_ind(end) - 1 + 2 * (1:length(dft_pos_ind))';
dfs_imag_ind = dfs_sym_ind(end) + 2 * (1:length(dft_pos_ind))';

% construct the transformation matrix from dft to dfs (the real fourier
% series)
dfs_matrix = dft2dfs_matrix(dft_sym_ind, dft_pos_ind, dft_neg_ind, dfs_sym_ind, dfs_real_ind, dfs_imag_ind);

% create vectorized desired correlation output
yf_vec = single(yf(:));

yf_vec_perm = permute(yf_vec(dft_nonneg_ind),[3 2 1]);

if params.use_reg_window
    % create weight window
    ref_window_power = params.reg_window_power;
    
    % normalization factor
    reg_scale = 0.5 * base_target_sz/featureRatio;
    
    % construct grid
    wrg = -(use_sz(1)-1)/2:(use_sz(1)-1)/2;
    wcg = -(use_sz(2)-1)/2:(use_sz(2)-1)/2;
    [wrs, wcs] = ndgrid(wrg, wcg);
    
    % construct the regukarization window
    reg_window = (params.reg_window_edge - params.reg_window_min) * (abs(wrs/reg_scale(1)).^ref_window_power + abs(wcs/reg_scale(2)).^ref_window_power) + params.reg_window_min;
    
    % compute the DFT and enforce sparsity
    reg_window_dft = fft2(reg_window) / prod(use_sz);
    reg_window_dft_sep = cat(3, real(reg_window_dft), imag(reg_window_dft));
    reg_window_dft_sep(abs(reg_window_dft_sep) < params.reg_sparsity_threshold * max(abs(reg_window_dft_sep(:)))) = 0;
    reg_window_dft = reg_window_dft_sep(:,:,1) + 1i*reg_window_dft_sep(:,:,2);
    
    % do the inverse transform, correct window minimum
    reg_window_sparse = real(ifft2(reg_window_dft));
    reg_window_dft(1,1) = reg_window_dft(1,1) - support_sz * min(reg_window_sparse(:)) + params.reg_window_min;
    
    % construct the regularizsation matrix
    regW = cconvmtx2(reg_window_dft);
    
    regW_dfs = real(dfs_matrix * regW * dfs_matrix');
    
    WW_block = regW_dfs' * regW_dfs;
    
    % If the filter size is small enough, remove small values in WW_block.
    % It takes too long time otherwise.
    if support_sz <= 120^2
        WW_block(0<abs(WW_block) & abs(WW_block)<0.00001) = 0;
    end
else
    % else use a scaled identity matrix
    WW_block = lambda * speye(support_sz);
    params.reg_window_min = sqrt(lambda);
end

% create block diagonal regularization matrix
WW = eval(['blkdiag(WW_block' repmat(',WW_block', 1, feature_dim-1) ');']);

WW_L = tril(WW);
WW_U = triu(WW, 1);

if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    
    scaleFactors = scale_step .^ scale_exp;
    
    %force reasonable scale changes
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

num_sym_coef = length(dft_sym_ind);

% create indexing vectors

% first create the indices for the symmetrix (real) part of the spectrum
index_i_sym = zeros(2*feature_dim, length(dft_sym_ind), feature_dim);
index_j_sym = zeros(size(index_i_sym));

index_i_sym_re = repmat(bsxfun(@plus, support_sz*(0:feature_dim-1)', 1:length(dft_sym_ind)), [1 1 feature_dim]); %index for the Real-Real part
index_i_sym(1:2:end, :, :) = index_i_sym_re;
index_i_sym(2:2:end, :, :) = NaN; % these will be zero

index_j_sym_re = permute(index_i_sym_re, [3 2 1]);
index_j_sym(1:2:end, :, :) = index_j_sym_re;
index_j_sym(2:2:end, :, :) = NaN; % these will be zero

% create the indices for the remaining part
index_i = zeros(2*feature_dim, 2*length(dft_pos_ind), feature_dim);
index_j = zeros(size(index_i));

index_i_re = repmat(bsxfun(@plus, support_sz*(0:feature_dim-1)', (length(dft_sym_ind)+1:2:support_sz)), [1 1 feature_dim]); %index for the Real-Real part
index_i(1:2:end, 1:2:end, :) = index_i_re;
index_i(2:2:end, 1:2:end, :) = index_i_re + 1;
index_i(1:2:end, 2:2:end, :) = index_i_re;
index_i(2:2:end, 2:2:end, :) = index_i_re + 1;

index_j_re = permute(index_i_re, [3 2 1]);
index_j(1:2:end, 1:2:end, :) = index_j_re;
index_j(2:2:end, 1:2:end, :) = index_j_re;
index_j(1:2:end, 2:2:end, :) = index_j_re + 1;
index_j(2:2:end, 2:2:end, :) = index_j_re + 1;

% concatenate the results
index_i = cat(2, index_i_sym, index_i);
index_j = cat(2, index_j_sym, index_j);

index_i = index_i(:);
index_j = index_j(:);

% the imaginary part of the autocorrelations (along the diagonal) will be zero
zero_ind = (index_i == index_j-1) | (index_i == index_j+1);
index_i(zero_ind) = NaN;
index_j(zero_ind) = NaN;

% indexing masks for upper and lower triangular part
data_L_mask = index_i >= index_j;
data_U_mask = index_i < index_j;

data_L_i = index_i(data_L_mask);
data_L_j = index_j(data_L_mask);
data_U_i = index_i(data_U_mask);
data_U_j = index_j(data_U_mask);

% extract the linear indeces from the data matrix and regularization matrix
WW_L_ind = find(WW_L);
% WW_U_ind = find(WW_U);
data_L_ind = sub2ind(size(WW_L), data_L_i, data_L_j);
% data_U_ind = sub2ind(size(WW_U), data_U_i, data_U_j);

% compute the linear indeces of the non-zeros in the matrix
[L_ind, ~, data_WW_in_L_index] = unique([data_L_ind; WW_L_ind]);
% [U_ind, ~, data_WW_in_U_index] = unique([data_U_ind; WW_U_ind]);

% compute the corresponding indices for the values in the data and reg
% matrix
data_in_L_index = uint32(data_WW_in_L_index(1:length(data_L_ind)));
% data_in_U_index = data_WW_in_U_index(1:length(data_U_ind));
WW_in_L_index = data_WW_in_L_index(length(data_L_ind)+1:end);
% WW_in_U_index = data_WW_in_U_index(length(data_U_ind)+1:end);

% create the arrays of values in the regularization matrix
nnz_L = length(L_ind);
% nnz_U = length(U_ind);
WW_L_vec = zeros(nnz_L, 1, 'single');
% WW_U_vec = zeros(nnz_U, 1);
WW_L_vec(WW_in_L_index) = full(WW_L(WW_L_ind));
% WW_U_vec(WW_in_U_index) = WW_U(WW_U_ind);

% precompute the data part of the regularization matrix
WW_L_vec_data = WW_L_vec(data_in_L_index);

% initialize the content vectors for the matrices
L_vec = WW_L_vec;
% U_vec = WW_U_vec;

% preallocate the matrices
mat_size = feature_dim * support_sz;
[L_i, L_j] = ind2sub(size(WW_L), L_ind);
% [U_i, U_j] = ind2sub(size(WW_U), U_ind);
AL = sparse(L_i, L_j, ones(nnz_L,1), mat_size, mat_size);
% AU = sparse(U_i, U_j, ones(nnz_U,1), mat_size, mat_size);
AU_data = sparse(data_U_i, data_U_j, ones(length(data_U_i),1), mat_size, mat_size);

if interpolate_response >= 3
    % Pre-computes the grid that is used for gradient ascent.
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    gradient_ascent_iterations = params.gradient_ascent_iterations;
end

rect_position = zeros(num_frames, 4);

disp_layer = 2;

time = 0;

xxlf_sep = zeros(2*feature_dim, length(dft_sym_ind) + 2 * length(dft_pos_ind), feature_dim, 'single');

multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales, 'uint8');

% initialization
prior_weights = [];
sample_weights = [];
latest_ind = [];
sample_frame = nan(params.nSamples);
dft_part_sz = cumsum([length(dft_sym_ind) length(dft_pos_ind) length(dft_neg_ind)]);
samplesf = 1i*zeros(params.nSamples,feature_dim,dft_part_sz(2),'single');
hf_rhs_reshaped = 1i*zeros(support_sz, feature_dim,'single');


for frame = 1:num_frames,
    %load image
    im = imread(s_frames{frame});
    if size(im,3) > 1 && colorImage == false
        im = im(:,:,1);
    end

    tic();
    
    if frame > 1
        old_pos = inf(size(pos));
        iter = 1;
        
        %translation search
        while iter <= refinement_iterations && any(old_pos ~= pos)
            % Get multi-resolution image
            for scale_ind = 1:nScales
                multires_pixel_template(:,:,:,scale_ind) = ...
                    get_pixels(im, pos, round(sz*currentScaleFactor*scaleFactors(scale_ind)), sz);
            end
            
            xt = bsxfun(@times,get_features(multires_pixel_template,features,global_feat_params),cos_window);
            
            xtf = fft2(xt);
            
            responsef = permute(sum(bsxfun(@times, hf, xtf), 3), [1 2 4 3]);
            
            if interpolate_response == 2
                % use dynamic interp size, resize to pixel coordinates
                interp_sz = floor(size(y) * featureRatio * currentScaleFactor);
            end
            responsef_padded = resizeDFT2(responsef, interp_sz);
            
            % response
            response = ifft2(responsef_padded, 'symmetric');
            
            % find maximum
            if interpolate_response >= 3
                [disp_row, disp_col, sind] = resp_newton(response, responsef_padded, gradient_ascent_iterations, ky, kx, use_sz);
            else
                [row, col, sind] = ind2sub(size(response), find(response == max(response(:)), 1));
                disp_row = mod(row - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
                disp_col = mod(col - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
            end
            
            % calculate translation
            switch interpolate_response
                case 0
                    translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
                case 1
                    translation_vec = round([disp_row, disp_col] * currentScaleFactor * scaleFactors(sind));
                case 2
                    translation_vec = round([disp_row, disp_col] * scaleFactors(sind));
                case {3, 4}
                    translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
            end
            
            % set the scale
            currentScaleFactor = currentScaleFactor * scaleFactors(sind);
            % adjust to make sure we are not to large or to small
            if currentScaleFactor < min_scale_factor
                currentScaleFactor = min_scale_factor;
            elseif currentScaleFactor > max_scale_factor
                currentScaleFactor = max_scale_factor;
            end
            
            % update position
            old_pos = pos;
            pos = pos + translation_vec;
            
            iter = iter + 1;
        end
        
    end
    
    % Update the prior weights
    [prior_weights, replace_ind] = update_prior_weights(prior_weights, sample_weights, latest_ind, frame, params);
    latest_ind = replace_ind;
    sample_frame(replace_ind) = frame;
    
    % Initialize the weight for the new sample
    if frame == 1
        sample_weights = prior_weights;
    else
       % ensure that the new sample always get its current prior weight
       new_sample_weight = learning_rate;
       sample_weights = sample_weights * (1 - new_sample_weight) / (1 - sample_weights(replace_ind));
       sample_weights(replace_ind) = new_sample_weight;
       sample_weights = sample_weights / sum(sample_weights);
    end;
    
    % Extract training sample
    pixels = get_pixels(im,pos,round(sz*currentScaleFactor),sz);
    xl = bsxfun(@times,get_features(pixels,features,global_feat_params),cos_window);
    
    xlf = fft2(xl);
    xlf_reshaped = reshape(xlf, [support_sz, feature_dim]);
    xlf_nonneg = xlf_reshaped(dft_nonneg_ind,:);
    
    % Store new sample
    samplesf(replace_ind,:,:) = permute(xlf_nonneg,[3,2,1]);
    
    for acs_iter = 1:params.num_acs_iter
        % Update sample weights (alpha in the paper)
        if ~params.sample_weight_after_filter
            if frame > params.sample_burnin && params.sample_reg > 0
                sample_loss = compute_sample_loss(hf_nonneg, samplesf, yf_vec_perm, dft_part_sz);
                sample_weights = update_weights(sample_loss, prior_weights, frame, params);
            else
                sample_weights = prior_weights;
            end;
        end
        
        % Data matrix of teh training samples
        sample_matrix = mtimesx(samplesf,'c',bsxfun(@times,sample_weights,samplesf),'speed');
        
        % compute autocorrelation
        xxlf_sym = permute(sample_matrix(:,:,1:dft_part_sz(1)),[1 3 2]);
        xxlf_pos = permute(sample_matrix(:,:,dft_part_sz(1)+1:dft_part_sz(2)),[1 3 2]);
        xxlf_pos_real = real(xxlf_pos);
        
        % partition the real and imaginary parts
        xxlf_sep(1:2:end, 1:num_sym_coef, :) = real(xxlf_sym);
        xxlf_sep(1:2:end, num_sym_coef+1:2:end, :) = xxlf_pos_real;
        xxlf_sep(2:2:end, num_sym_coef+1:2:end, :) = imag(xxlf_pos);
        xxlf_sep(1:2:end, num_sym_coef+2:2:end, :) = -imag(xxlf_pos);
        xxlf_sep(2:2:end, num_sym_coef+2:2:end, :) = xxlf_pos_real;
        
        % The data part of the left-hand side
        hf_autocorr = xxlf_sep(:);
        
        % The right-hand side
        weighted_samples = permute(mtimesx(sample_weights',samplesf,'speed'),[3,2,1]);
        hf_rhs_reshaped(dft_sym_ind,:) = weighted_samples(1:dft_part_sz(1),:);
        hf_rhs_reshaped(dft_pos_ind,:) = weighted_samples(dft_part_sz(1)+1:end,:);
        hf_rhs_reshaped(dft_neg_ind,:) = conj(weighted_samples(dft_part_sz(1)+1:end,:));
        hf_rhs_reshaped = bsxfun(@times,yf_vec, conj(hf_rhs_reshaped));
        hf_rhs_reshaped_real = real(dfs_matrix * double(hf_rhs_reshaped));
        hf_rhs = hf_rhs_reshaped_real(:);
        
        % Initialization of the filter
        if frame == 1 && acs_iter == 1
            hf_init_autocorr = double(sum(xlf_reshaped .* conj(xlf_reshaped), 2));
            switch params.init_strategy
                case 'const_reg'       % exact solution for constant regularization
                    hf_init = bsxfun(@rdivide, xyf_corr, hf_init_autocorr + params.reg_window_min^2);
                    hf_init = real(dfs_matrix * hf_init);
                    hf_vec = hf_init(:);
                case 'indep'           % independent filters for each feature
                    A_init = real(dfs_matrix * spdiags(hf_init_autocorr, 0, support_sz, support_sz) * dfs_matrix') + feature_dim * WW_block;
                    b_init = reshape(hf_rhs, support_sz, feature_dim);
                    hf_init = A_init \ b_init;
                    hf_vec = hf_init(:);
            end
            
            num_GS_iter = params.num_init_GS_iter;
            
            if debug == 1
                hf_vec_old = hf_vec;
                hf_difference = zeros(num_frames * num_GS_iter, 1);
                
                % reconstruct the filter
                hf = reshape(single(dfs_matrix' * reshape(hf_vec, [support_sz, feature_dim])), [use_sz, feature_dim]);
                joint_loss = []; hf_grad_norm = []; weights_grad_norm = [];
            end
        else
            num_GS_iter = params.num_GS_iter;
        end
        
        if debug
            [joint_loss(end+1), hf_grad_norm(end+1), weights_grad_norm(end+1)] = eval_optimality(hf, sample_weights, prior_weights, samplesf, regW, dft_nonneg_ind, dft_part_sz, yf_vec_perm, frame, params);
        end
        
        % add the autocorrelation to the matrix vectors with the regularization
        L_vec(data_in_L_index) = hf_autocorr(data_L_mask) + WW_L_vec_data;
        
        % update the matrices with the new non-zeros
        AL = setnonzeros(AL, double(L_vec));
        AU_data = setnonzeros(AU_data, double(hf_autocorr(data_U_mask)));
        
        % do Gausss-Seidel
        for iter = 1:num_GS_iter
            hf_vec = AL \ (hf_rhs - AU_data * hf_vec - WW_U * hf_vec);
            
            if debug
                hf_difference((frame - 1)*num_GS_iter + iter) = norm(hf_vec - hf_vec_old);
                hf_vec_old = hf_vec;
            end
        end
        
        % reconstruct the filter
        hf = reshape(single(dfs_matrix' * reshape(hf_vec, [support_sz, feature_dim])), [use_sz, feature_dim]);
        
        if debug
            [joint_loss(end+1), hf_grad_norm(end+1), weights_grad_norm(end+1)] = eval_optimality(hf, sample_weights, prior_weights, samplesf, regW, dft_nonneg_ind, dft_part_sz, yf_vec_perm, frame, params);
        end
        
        % compute reshaped filter for using in the weight determination
        hf_reshaped = reshape(hf, [support_sz, feature_dim]);
        hf_nonneg = hf_reshaped(dft_nonneg_ind,:);
        
        % Update sample weights (alpha in the paper)
        if params.sample_weight_after_filter
            if frame > params.sample_burnin && params.sample_reg > 0
                sample_loss = compute_sample_loss(hf_nonneg, samplesf, yf_vec_perm, dft_part_sz);
                sample_weights = update_weights(sample_loss, prior_weights, frame, params);
            else
                sample_weights = prior_weights;
            end;
        end
    end
    
    if debug
        [joint_loss(end+1), hf_grad_norm(end+1), weights_grad_norm(end+1)] = eval_optimality(hf, sample_weights, prior_weights, samplesf, regW, dft_nonneg_ind, dft_part_sz, yf_vec_perm, frame, params);
    end
    
    % debug visualization
    if debug
        history_len = 100;
        fontSizeAxis = 12;
        
        % Plot joint cost
        figure(44);
        si = max(1, length(joint_loss) - history_len);
        plot(joint_loss(si:end), 'b-x', 'linewidth',1.5); title('Joint Loss');
        axis([1 length(joint_loss(si:end)) min(joint_loss(si:end)) max(joint_loss(si:end))]);
        set(gca, 'fontSize', fontSizeAxis);
        
        % Plot optimality measure
        figure(46);
        total_gradient_norm = sqrt(hf_grad_norm.^2 + weights_grad_norm.^2);
        si = max(1, length(total_gradient_norm) - history_len);
        semilogy(total_gradient_norm(si:end), 'b.-', 'linewidth',1, 'markersize', 15); title('KKT Optimality Measure');
        hold on;
        plot(hf_grad_norm(si:end), 'g--', 'linewidth',2);
        plot(weights_grad_norm(si:end), 'r--', 'linewidth',2); hold off;
        axis([1 length(total_gradient_norm(si:end)) 0.1*min(total_gradient_norm(si:end)) 1.1*max(total_gradient_norm(si:end))]);
        ylabel('Lagrangian gradient magnitude')
        legend({'Total gradient', 'Filter parameters gradient', 'Weight parameters gradient'}, 'location', 'southeast')
        set(gca, 'fontSize', fontSizeAxis);
        
        % Sample weights
        figure(71);
        plot_frames = nan(num_frames,1); plot_sample = nan(num_frames,1); plot_prior = nan(num_frames,1);
        max_ind = min(frame, params.nSamples);
        [sorted_frames, ind] = sort(sample_frame(1:max_ind));
        plot_frames(sorted_frames) = sorted_frames;
        plot_sample(sorted_frames) = sample_weights(ind);
        plot_prior(sorted_frames) = prior_weights(ind);
        plot(plot_frames, plot_prior, 'or-', 'linewidth',1.5, 'markersize', 3);
        hold on;
        plot(plot_frames, plot_sample, 'xb-', 'linewidth',1.5, 'markersize', 3);
        title('Sample and prior weights');
        xlabel('Frame number');
        legend({'Prior weights', 'Sample weights'}, 'location', 'southeast')
        hold off;
        axis([1 num_frames 0 1.1*max([sample_weights; prior_weights])]);
        set(gca, 'fontSize', fontSizeAxis);
        
        % Filter
%         figure(20)
%         subplot_cols = 1;%ceil(sqrt(feature_dim));
%         subplot_rows = 1;%ceil(feature_dim/subplot_cols);
%         for disp_layer = 1:1;%feature_dim
%             subplot(subplot_rows,subplot_cols,disp_layer);
%             imagesc(ifft2(conj(hf(:,:,disp_layer)), 'symmetric')); 
%             colorbar;
%             axis image;
%         end
%         figure(99);plot(hf_difference);axis([1, num_GS_iter * frame, 0, max(hf_difference)]);
    end
    
    target_sz = floor(base_target_sz * currentScaleFactor);
    
    %save position and calculate FPS
    rect_position(frame,:) = [pos([2,1]) - floor(target_sz([2,1])/2), target_sz([2,1])];
    
    time = time + toc();
    
    %visualization
    if params.visualization == 1
        rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        im_to_show = double(im)/255;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end
        if frame == 1,  %first frame, create GUI
            fig_handle = figure('Name', 'Tracking');
%             set(fig_handle, 'Position', [100, 100, size(im,2), size(im,1)]);
            imagesc(im_to_show);
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            hold off;
            axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])

        else
            try  %subsequent frames, update GUI
                resp_sz = round(sz*currentScaleFactor*scaleFactors(scale_ind));
                xs = floor(old_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
                ys = floor(old_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);
                sc_ind = floor((nScales - 1)/2) + 1;
                
                figure(fig_handle);
%                 set(fig_handle, 'Position', [100, 100, 100+size(im,2), 100+size(im,1)]);
                imagesc(im_to_show);
                hold on;
                resp_handle = imagesc(xs, ys, fftshift(response(:,:,sc_ind))); colormap hsv;
                alpha(resp_handle, 0.5);
                rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
                text(10, 10, int2str(frame), 'color', [0 1 1]);
                axis off;axis image;
                hold off;

%                 axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
                
            catch
                return
            end
        end
        
        drawnow
    end
end

% close(writer);

fps = numel(s_frames) / time;

disp(['fps: ' num2str(fps)])

results.type = 'rect';
results.res = rect_position;%each row is a rectangle
results.fps = fps;