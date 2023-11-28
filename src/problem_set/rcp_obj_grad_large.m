function [f, g_h] = rcp_obj_grad_large(R, X,flag)

% compute function value and its gradient when f = <C,R'*R>, and C = -X*X';



if nargin<3; flag=1; end
% 


% function value
if(flag==1)
    RC = -(R * X)*X';
     f = sum(sum(R .* RC));
else
    f = 0;
end
if nargout > 1
    % gradient handle
    g_h = @(V) -(V * X)*X';
end

end