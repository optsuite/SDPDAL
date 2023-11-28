function [f, g_h, g] = spca_obj_grad(R, C,flag)
if nargin<3; flag=1; end
% compute RC


% function value
if(flag==1)
    RC = (R * C);
f = sum(sum(R .* RC));
else
    f = 0;
end
if nargout > 1
    % gradient handle
    g_h = @(V) (V * C);
end

% real g
if nargout > 2
    g = (R * C);
end

end