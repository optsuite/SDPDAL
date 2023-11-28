function fobj = mc_factory(type, C, alpha, lambda)
fobj = struct();
fobj.data = {};
if strcmp(type, 'T')
    if mod(alpha, 1) ~= 0 || alpha <= 1
        error('alpha must be integer > 1');
    end
    fobj.obj_grad = @(R) fun_entropy_T_int(R, C, alpha, lambda);
    fobj.hess = @(R, U, data) hess_entropy_T_int(R, U, data, alpha, lambda);
    fobj.fun_extra = @(data) prepare_cache(data, alpha);
elseif strcmp(type, 'R')
    if mod(alpha, 1) ~= 0 || alpha <= 1
        error('alpha must be integer > 1');
    end
    fobj.obj_grad = @(R) fun_entropy_R_int(R, C, alpha, lambda);
    fobj.hess = @(R, U, data) hess_entropy_R_int(R, U, data, alpha, lambda);
    fobj.fun_extra = @(data) prepare_cache(data, alpha);
elseif strcmp(type, 'N')
    fobj.obj_grad = @(R) fun_plain(R, C);
    fobj.hess = @(~, ~, ~) 0;
else
    error(['unknown type ' type]);
end
end

function [f, g_h, g] = fun_plain(R, C)
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
end

function [f, g_h, g] = fun_entropy_T_int(R, C, alpha, lambda)
    % suppose alpha is an integer
    RRt = R * R'; % should be small
    RRtam2 = RRt ^ (alpha - 2); % matrix power alpha - 2
    RRtam1 = RRt * RRtam2;      % matrix power alpha - 1
     
    trX = trace(RRt);
    trXa = sum(sum(RRt .* RRtam1));
    
    % obj
    g1 = R * C;
    f1 = sum(sum(g1 .* R));
    f2 = trXa / (trX ^ alpha);
    f = f1 + lambda * (f2 - 1) / (1 - alpha);
    
    % gradient handle
    if nargout > 1
        RXam2 = RRtam2 * R;
        g2_h = @(V) (V * R') * RXam2 - trXa / trX * V;
        g_h = @(V) V * C + lambda * alpha / (1 - alpha) / (trX ^ alpha) * g2_h(V);
    end
    
    % true gradient
    if nargout > 2
        g2 = RRt * RXam2 - trXa / trX * R;
        g = g1 + lambda * alpha / (1 - alpha) / (trX ^ alpha) * g2;
    end
end

function h = hess_entropy_T_int(R, U, data, alpha, lambda)
    % load from cache
    RRt = data.RRt;
    RXam1 = data.RXam1;
    trX = data.trX;
    trXa = data.trXa;
    
    % comput RU'
    RUt = R * U';
    
    % second part
    h2 = -2 * alpha * sum(sum(U .* R)) / trX ^ (alpha + 1) * RXam1 + ...
        2 * ((alpha + 1) * trXa - alpha * sum(sum(U .* RXam1))) / trX ^ (alpha + 2) * R;
    
    % third part

    h3 = RRt * U + RUt * R;
    h3_plus = h3;
    for i = 2:(alpha - 1)
        h3_plus = RRt * h3_plus;
        h3 = h3_plus + (h3 * R') * R;
    end
    
    % combine these terms
    h = lambda * alpha / (1 - alpha) * (h2 + h3 / trX ^ alpha);
end

function [f, g_h, g] = fun_entropy_R_int(R, C, alpha, lambda)
    % suppose alpha is an integer
    RRt = R * R'; % should be small
    RRtam2 = RRt ^ (alpha - 2); % matrix power alpha - 2
    RRtam1 = RRt * RRtam2;      % matrix power alpha - 1
    
    trX = trace(RRt);
    trXa = sum(sum(RRt .* RRtam1));
    
    % obj
    g1 = R * C;
    f1 = sum(sum(g1 .* R));
    f2 = log(trXa / (trX ^ alpha));
    f = f1 + lambda * f2 / (1 - alpha);
    
    % gradient handle
    if nargout > 1
        RXam2 = RRtam2 * R;
        g2_h = @(V) (V * R') * RXam2 / trXa - V / trX;
        g_h = @(V) V * C + lambda * alpha / (1 - alpha) * g2_h(V);
    end
    
    % true gradient
    if nargout > 2
        g2 = RRt * RXam2 / trXa - R / trX;
        g = g1 + lambda * alpha / (1 - alpha) * g2;
    end
end

function h = hess_entropy_R_int(R, U, data, alpha, lambda)
    % load from cache
    RRt = data.RRt;
    RXam1 = data.RXam1;
    trX = data.trX;
    trXa = data.trXa;
    
    % comput RU'
    RUt = R * U';
    
    % second part
    h2 = -2 * alpha * sum(sum(U .* RXam1)) / trXa ^ 2 * RXam1 + ...
        2 * sum(sum(U .* R)) / trX ^ 2 * R;
    
    % third part

    h3 = RRt * U + RUt * R;
    h3_plus = h3;
    for i = 2:(alpha - 1)
        h3_plus = RRt * h3_plus;
        h3 = h3_plus + (h3 * R') * R;
    end
    
    % combine these terms
    h = lambda * alpha / (1 - alpha) * (h2 + h3 / trXa);
end

function data = prepare_cache(data, alpha)
    R = data.XP;
    
    data.RRt = R * R';
    RRtam1 = data.RRt ^ (alpha - 1);
    data.RXam1 = RRtam1 * R;
    data.trX = trace(data.RRt);
    data.trXa = sum(sum(data.RRt .* RRtam1));
end
