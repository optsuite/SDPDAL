function [f, g_h, g] = nomad_obj_grad(R, Y)
% rcp function value
% f(R'R) = <YY', R'R>

% compute RY
RY = R * Y;


% function value
f = -norm(RY, 'fro') ^ 2;

if nargout > 1
    % gradient handle
    g_h = @(V) -(V * Y)*Y';
end

if nargout > 2
    % real gradient (note: intermediate results can be used)
    g = -RY * Y';
end