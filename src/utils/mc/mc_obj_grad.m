function [f, g_h, g] = mc_obj_grad(R, C)
% compute RC
RC = R * C;

% function value
f = sum(sum(R .* RC));

if nargout > 1
    % gradient handle
    g_h = @(V) (V * C);
end

if nargout > 2
    % gradient
    g = RC;
end