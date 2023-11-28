function ret = EPSDP_MC(C, opts)

n = size(C, 1);

% default parameters
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'p0'); opts.p0 = 0; end
if ~isfield(opts, 'entropy'); opts.entropy = 'tsallis'; end
if ~isfield(opts, 'alpha'); opts.alpha = 2; end
if ~isfield(opts, 'lambda0'); opts.lambda0 = 10; end
if ~isfield(opts, 'verbose'); opts.verbose = 1; end
if ~isfield(opts, 'maxit'); opts.maxit = 30; end
if ~isfield(opts, 'tol'); opts.tol = 1e-3; end
if ~isfield(opts, 'rho'); opts.rho = 1; end
if ~isfield(opts, 'second_order'); opts.second_order = 1; end

if opts.p0 == 0
    p0 = min(ceil(sqrt(2 * n)), 200);
else
    p0 = opts.p0;
end

% initialization
x0 = randn(p0,n); nrmx0 = dot(x0,x0,1);
x0 = bsxfun(@rdivide, x0, sqrt(nrmx0));
x = x0;
p = p0;
M = obliquefactory(p,n);

% regularization parameters
alpha = opts.alpha;
lambda = opts.lambda0;

% use gradient method for inital iterations
do_NEWT = 0;

% default options for arnt
optARNT = default_ARNT_options(opts);

% string format (for logging)
str_header = 'iter       obj          cut        nrmG    rank    lambda   ratio-1  time\n';
str_header_l = '-----------------------------------------------------------------------------------\n';
str_data = '%4d %12.6e % 12.6e %6.2e %5d   %6.2e %6.2e  %.2f';

tt = tic;
% main loop
for iter = 1:opts.maxit
    % perpare iterations for arnt
    if strcmp(opts.entropy, 'tsallis')
        func_handle = @fun_entropy_T;
        grad_handle = @(x) grad_entropy_T(x, C, alpha, lambda);
        optARNT.fun_extra = @(data) prepare_cache(data, alpha);
        optARNT.hess = @(R, U, data) hess_entropy_T(R, U, data, C, alpha, lambda);
    elseif strcmp(opts.entropy, 'renyi')
        func_handle = @fun_entropy_R;
        grad_handle = @(x) grad_entropy_R(x, C, alpha, lambda);
        optARNT.fun_extra = @(data) prepare_cache(data, alpha);
        optARNT.hess = @(R, U, data) hess_entropy_R(R, U, data, C, alpha, lambda);
    elseif strcmp(opts.entropy, 'none') || iter == 0 % for plain SDP, or the initialization for entropy
        func_handle = @fun;
        grad_handle = @(x) grad(x, C);
        optARNT.hess = @(X, U, data) hess(X, U, data, C);
        optARNT.fun_extra = @(data) fun_extra(data);
        optARNT.maxit = 500;
    else
        error('unknown entropy type %s', opts.entropy);
    end
    optARNT.grad = grad_handle;
    
    if strcmp(opts.entropy, 'none')
        [x, ~, out] = arnt(x0, @fun, M, optARNT, C);
        time = out.time;
    elseif do_NEWT == 1
        % call arnt
        [x, ~, out] = arnt(x, func_handle, M, optARNT, C, alpha, lambda);
        time = out.time;
    else
        % call RGBB (do not invoke second order method)
        t = tic;
        [x, ~, out] = RGBB(x, func_handle, M, optARNT.opts_init, C, alpha, lambda);
        time = toc(t);
    end
    
    % now check the rank of x
    [U, D, V] = svd(x, 'econ');
    d = diag(D);
    rank_x = sum(d / d(1) > 1e-5);
    new_p = max(3, min(p0, ceil(1.1 * rank_x)));
    if new_p < p
        x = U(:, 1:new_p) * diag(d(1:new_p)) * V(:, 1:new_p)';
    end
    
    % compute ratio
    ratio = sum(d) / d(1);
    
    if opts.verbose
        if iter == 1
            fprintf(str_header);
            fprintf(str_header_l);
        end
        cut = fun(x, C);
        fprintf(str_data, iter, out.fval, cut, out.nrmG, p, lambda, ...
            ratio - 1, time);
        if new_p < p
            fprintf('  p:[%d -> %d]', p, new_p);
        end
        fprintf('\n');
    end
    
    % termination rule
    if ratio - 1 < opts.tol && opts.second_order && do_NEWT == 1 || strcmp(opts.entropy, 'none')
        break;
    end
    
    % switch to second_order algorithm when ratio - 1 is small
    if ratio - 1 < 1e2 * opts.tol
        if opts.second_order == 1 && mod(alpha, 1) == 0
            do_NEWT = 1;
        end
    end
    
    % update lambda and p
    lambda = ratio ^ opts.rho * lambda;
    p = new_p;
end

ret = struct();
ret.time = toc(tt);
ret.R = x;
ret.iter = iter;
ret.obj = cut;
ret.pinf = norm(dot(x,x,1) - 1, inf);
ret.nrmG = out.nrmG;
ret.p = p;

end

function options = default_ARNT_options(opts)

options = struct();
options.record = 1;
options.gtol = 1e-6;
options.xtol = 1e-5;
options.ftol = 1e-14;
%M = obliquefactory(p,n);
%optARNT.hess = @(X, U, data) hess(X, U, data, C);
%optARNT.grad = @(X) grad(X, C);
%options.grad = @(X) grad_entropy_T(X, C, alpha, lambda);
options.fun_TR = @fun_sub;
options.fun_extra = @fun_extra;
options.opts_init.record = opts.verbose > 1;
options.solver_init = @RGBB;
options.opts_init.tau   = 1e-6;
options.opts_init.maxit = 3000;
options.opts_init.gtol  = options.gtol*1e4;
options.opts_init.xtol  = options.xtol*1e2;
options.opts_init.ftol  = options.ftol*1e2;
options.opts_sub.record = 0;
options.solver_sub  = @RNewton;
options.opts_sub.record = 0;
options.opts_sub.tau    = 1e-6;
options.opts_sub.maxitvec  = [100,150,200,300,1000];
options.opts_sub.gtol   = options.gtol*1e0;
options.opts_sub.xtol   = options.xtol*1e0;
options.opts_sub.ftol   = options.ftol*1e0;
options.fun_TR = @fun_sub;
options.maxit = 1000; 
options.tau = 10;
options.usenumstab = 0;

end


function [f,g] = fun(X,C)
    % max Tr(C*Z), s.t., Z_ii = 1, Z psd
    % low rank model:
    % Z = X'*X, X = [X_1, ..., X_n], X is a p by n matrix
    % max Tr(C*X'*X), s.t., ||X_i|| = 1,
    g = 2*(X*C);
    f = sum(dot(g,X))/2;
end

function data = fun_extra(data)
    %XP = data.XP;
end

function g = grad(X, C)
    g = 2*(X*C);
end

function h = hess(~, U, data, C)
    h = 2*U*C + data.sigma*U;
end

%     function [f,g] = fun_sub(X,data)
%         tempc = 2*(X*C);
%         UmX = X - data.XP;
%         nrm = norm(UmX,'fro');
%         f = .5*sum(sum(X.*tempc))+1/3*data.TrRho*nrm^3;
%         g = tempc + data.TrRho*nrm*UmX;
%         
%         %f = .5*sum(sum(X.*tempc))+1/2*data.TrRho*nrm^2;
%         %g = tempc + data.TrRho*UmX;        
%     end

function [f, g] = fun_entropy_T(R, C, alpha, lambda)
    % X = R'R
    % R = UDV' with V'V = I
    if mod(alpha, 1) == 0
        [f, g] = fun_entropy_T_int(R, C, alpha, lambda);
        return
    end
    [U, D, V] = svd(R, 'econ');
    d = diag(D);
    
    trX = norm(d, 2) ^ 2;
    trXa = norm(d, 2 * alpha) ^ (2 * alpha);
    
    g1 = R * C;
    f1 = sum(dot(g1,R));
    f2 = trXa / (trX ^ alpha);
    f = f1 + lambda * (f2 - 1) / (1 - alpha);
    
    g2 = U * diag(d .^ (2 * alpha - 1)) * V' - trXa / trX * R;
    g = 2 * g1 + 2 * lambda * alpha / (1 - alpha) / (trX ^ alpha) * g2;
end

function [f, g] = fun_entropy_T_int(R, C, alpha, lambda)
    % suppose alpha is an integer
    RRt = R * R'; % should be small
    RRtam1 = RRt ^ (alpha - 1); % matrix power alpha - 1
    
    trX = trace(RRt);
    trXa = sum(sum(RRt .* RRtam1));
    RXam1 = RRtam1 * R;
    
    g1 = R * C;
    f1 = sum(dot(g1,R));
    f2 = trXa / (trX ^ alpha);
    f = f1 + lambda * (f2 - 1) / (1 - alpha);
    
    g2 = RXam1 - trXa / trX * R;
    g = 2 * g1 + 2 * lambda * alpha / (1 - alpha) / (trX ^ alpha) * g2;
end

function g = grad_entropy_T(R, C, alpha, lambda)
    [~, g] = fun_entropy_T(R, C, alpha, lambda);
end

function h = hess_entropy_T(R, U, data, C, alpha, lambda)
    if mod(alpha, 1) == 0
       h = hess_entropy_T_int(R, U, data, C, alpha, lambda);
       return
    end
    error('Not implemented!');
end

function h = hess_entropy_T_int(R, U, data, C, alpha, lambda)
    % load from cache
    RRt = data.RRt;
    RXam1 = data.RXam1;
    RXam2 = data.RXam2;
    trX = data.trX;
    trXa = data.trXa;
    
    % comput RU'
    RUt = R * U';
    
    % first part
    coeff1 = 2 * lambda * alpha  / (1 - alpha) / trX ^ alpha;
    coeff2 = data.sigma - 2 * lambda * alpha  / (1 - alpha) * trXa / trX ^ (alpha + 1);
    h1 = 2*U*C + coeff1 * RUt' * RXam2 + coeff2 * U;
    
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
    h = h1 + 2 * lambda * alpha / (1 - alpha) * (h2 + h3 / trX ^ alpha);
end

function [f, g] = fun_entropy_R(R, C, alpha, lambda)
    % X = R'R
    % R = UDV' with V'V = I
    if mod(alpha, 1) == 0
       [f, g] = fun_entropy_R_int(R, C, alpha, lambda);
       return
    end
    [U, D, V] = svd(R, 'econ');
    d = diag(D);
    
    trX = norm(d, 2) ^ 2;
    trXa = norm(d, 2 * alpha) ^ (2 * alpha);
    
    g1 = R * C;
    f1 = sum(dot(g1,R));
    f2 = log(trXa / (trX ^ alpha));
    f = f1 + lambda * f2 / (1 - alpha);
    
    g2 = U * diag(d .^ (2 * alpha - 1)/ trXa) * V' - R / trX;
    g = 2 * g1 + 2 * lambda * alpha / (1 - alpha) * g2;
end

function [f, g] = fun_entropy_R_int(R, C, alpha, lambda)
    % suppose alpha is an integer
    RRt = R * R'; % should be small
    RRtam1 = RRt ^ (alpha - 1); % matrix power alpha - 1
    
    trX = trace(RRt);
    trXa = sum(sum(RRt .* RRtam1));
    RXam1 = RRtam1 * R;
    
    g1 = R * C;
    f1 = sum(dot(g1,R));
    f2 = log(trXa / (trX ^ alpha));
    f = f1 + lambda * f2 / (1 - alpha);
    
    g2 = RXam1 / trXa - R / trX;
    g = 2 * g1 + 2 * lambda * alpha / (1 - alpha) * g2;
end

function g = grad_entropy_R(R, C, alpha, lambda)
    [~, g] = fun_entropy_R(R, C, alpha, lambda);
end

function h = hess_entropy_R(R, U, data, C, alpha, lambda)
    if mod(alpha, 1) == 0
       h = hess_entropy_R_int(R, U, data, C, alpha, lambda);
       return
    end
    error('Not implemented!');
end

function h = hess_entropy_R_int(R, U, data, C, alpha, lambda)
    % load from cache
    RRt = data.RRt;
    RXam1 = data.RXam1;
    RXam2 = data.RXam2;
    trX = data.trX;
    trXa = data.trXa;
    
    % comput RU'
    RUt = R * U';
    
    % first part
    coeff1 = 2 * lambda * alpha  / (1 - alpha) / trXa;
    coeff2 = data.sigma - 2 * lambda * alpha  / (1 - alpha) / trX;
    h1 = 2*U*C + coeff1 * RUt' * RXam2 + coeff2 * U;
    
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
    h = h1 + 2 * lambda * alpha / (1 - alpha) * (h2 + h3 / trXa);
end

function data = prepare_cache(data, alpha)
    R = data.XP;
    
    if mod(alpha, 1) ~= 0
        return
    end
    data.RRt = R * R';
    RRtam2 = data.RRt ^ (alpha - 2);
    RRtam1 = data.RRt * RRtam2;
    data.RXam1 = RRtam1 * R;
    data.RXam2 = RRtam2 * R;
    data.trX = trace(data.RRt);
    data.trXa = sum(sum(data.RRt .* RRtam1));
end