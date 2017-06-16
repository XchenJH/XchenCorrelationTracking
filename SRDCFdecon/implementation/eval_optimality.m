function [joint_loss, hf_grad_norm, weights_grad_norm] = eval_optimality(hf, sample_weights, prior_weights, samplesf, regW, dft_nonneg_ind, dft_part_sz, yf_vec_perm, frame, params )

% Evaluate optimality

[sz1, sz2, feature_dim] = size(hf);
support_sz = sz1*sz2;

norm_factor = 1 / support_sz;

% compute number of existing samples
num_samples = min(frame, params.nSamples);

sample_weights = sample_weights(1:num_samples);
prior_weights = prior_weights(1:num_samples);
samplesf = samplesf(1:num_samples,:,:);

%% --------- Gradient of the filter (theta) ---------
% compute reshaped filter
hf_reshaped = reshape(hf, [support_sz, feature_dim]);
hf_nonneg = hf_reshaped(dft_nonneg_ind,:);

% data term
corr_samp = mtimesx(samplesf, permute(hf_nonneg,[2 3 1]), 'speed');
abf = bsxfun(@times,sample_weights,corr_samp);
babf = mtimesx(samplesf,'c',abf,'speed');
hf_data_nonneg = permute(babf, [3 1 2]);

% regularization term
W_hf = regW * double(hf_reshaped);
hf_reg = single(regW' * W_hf);
hf_reg_nonneg = hf_reg(dft_nonneg_ind,:);

% rhs
weighted_samples = mtimesx(sample_weights',samplesf,'speed');
rhs_nonneg = permute(bsxfun(@times, conj(weighted_samples), yf_vec_perm), [3 2 1]);

hf_grad_nonneg = 2*(hf_data_nonneg + hf_reg_nonneg - rhs_nonneg) * norm_factor;

% hf gradient norm
hf_grad_norm = sqrt(conjugate_scalar_product(hf_grad_nonneg, hf_grad_nonneg, dft_part_sz, feature_dim));


%% --------- Gradient of the weights (alpha) ---------
    
% Compute residuals fo each frame
corr_error = bsxfun(@minus,corr_samp,yf_vec_perm);
residuals = norm_factor * (2*real(sum(corr_error .* conj(corr_error),3)) - ...
    real(sum(corr_error(:,1,1:dft_part_sz(1)) .* conj(corr_error(:,1,1:dft_part_sz(1))), 3)));

% Compute the usual gradient for the weights (alpha)
weights_grad = double(residuals + 2/params.sample_reg * (sample_weights ./ prior_weights));

% Compute KKT gradient
weights_zero_ind = (sample_weights < 1e-6);
lambda_zero_ind = ~weights_zero_ind;

% Matrices for quadprog
H = double([num_samples, 2*ones(1,num_samples); zeros(num_samples,1), eye(num_samples)]);
f = -[sum(weights_grad); weights_grad];
A = [zeros(num_samples,1), -eye(num_samples)];
b = zeros(num_samples,1);
Aeq = diag(double(lambda_zero_ind));
Aeq = [zeros(sum(lambda_zero_ind),1), Aeq(lambda_zero_ind,:)];
beq = zeros(sum(lambda_zero_ind),1);

options.Display = 'off';

lambda = quadprog(H, f, A, b, Aeq, beq,[],[],[],options);

weights_kkt_grad = weights_grad - lambda(2:end) - lambda(1)*ones(num_samples,1);

weights_grad_norm = norm(weights_kkt_grad);

%% --------- Joint loss ---------

joint_loss = sample_weights'*residuals + ...
    1/params.sample_reg * sum(sample_weights.^2 ./ prior_weights) + ...
    conjugate_scalar_product(W_hf, W_hf, dft_part_sz, feature_dim) * norm_factor;
