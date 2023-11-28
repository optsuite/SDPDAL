function [f, g_h, g] = theta_obj_grad(R)

n = size(R, 2);
% for theta problem, C = -E
sumR = sum(R, 2);
f = -norm(sumR, 2) ^ 2 / n;

% g_h is returned as a function handle
if nargout > 1
    g_h = @(V) repmat(-sum(V, 2), [1, n]) / n;
end

% real g
if nargout > 2
    g = repmat(-sumR, [1, n]) / n;
end