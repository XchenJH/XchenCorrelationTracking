function [ product ] = conjugate_scalar_product( a, b, sym_indices, feature_dim )
%CONJUGATE_SCALAR_PRODUCT Scalar product for conjugate symmetric signals
% computes a'*b when a and b are complex conjugate symmetric

ind = bsxfun(@plus, (1:sym_indices(1))', sym_indices(2) * (0:feature_dim-1));
product = 2*real(a(:)' * b(:)) - real(a(ind(:))' * b(ind(:)));

end

