function [f, g_h, g] = ncm_obj_grad(R, H, G)
% ncm function value
% f(R'R) = 1/2||H.*(X-G)||_F^2



%HHG = (H.^2).*(R'*R - G);
RR = R'*R;
if(~isempty(H))
    HG = H.*(RR - G);
    HHG = H.*HG;
else
    HG = RR - G;
    HHG = HG;
end
%HHG = HHG - mu*(H.^2).*RR;
% function value
f = 0.5*norm(HG, 'fro')^2 ;%- 0.5*mu*norm(H.*RR,'fro')^2;



if nargout > 1
    % gradient handle
    g_h = @(V) V * HHG;
end

if nargout > 2
    % real gradient (note: intermediate results can be used)
    grad = R * HHG;
    g.RR = RR; g.grad = grad;
end



