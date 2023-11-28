function [R, M] = do_decrease_p(new_p, R0, normR, man_handle)
% obtain size
n = size(R0, 2);
p = new_p;

% perform svd
[~, D, V] = svd(R0, 'econ');
R = D(1:p, 1:p) * V(:, 1:p)';

% normalize
R = R/sqrt(sum(sum(R.^2)))*normR;

% generate new manifold
M = man_handle(p,n,normR);

end