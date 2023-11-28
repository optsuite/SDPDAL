function ret = EPSDP_MC_ALM(Corig, opts)

n = size(Corig, 1);

% default parameters
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'p0'); opts.p0 = 0; end
if ~isfield(opts, 'entropy'); opts.entropy = 'T'; end
if ~isfield(opts, 'alpha'); opts.alpha = 2; end
if ~isfield(opts, 'lambda0'); opts.lambda0 = 0.001; end
if ~isfield(opts, 'verbose'); opts.verbose = 1; end
if ~isfield(opts, 'maxit'); opts.maxit = 50; end
if ~isfield(opts, 'tol'); opts.tol = 1e-3; end
if ~isfield(opts, 'rho'); opts.rho = 1; end
if ~isfield(opts, 'cuts'); opts.cuts = []; end

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

rho = opts.rho;

% processing cuts
if isempty(opts.cuts)
    AIop = [];
    bI = [];
    tol_max = 1e-4;
    tol_min = 1e-6;
else
    AIop = mc_gen_AIop(opts.cuts);
    bI = ones(numel(opts.cuts), 1);
    tol_max = 1e-3;
    tol_min = 1e-5;
end

% regularization parameters
alpha = opts.alpha;
lambda = opts.lambda0;

% default options for arnt
optSDPDAL = default_SDPDAL_options(opts);
if ~isempty(opts.cuts)
    optSDPDAL.comp_eigS = 2;
end

% scaling
normC = max(1, norm(Corig, 'fro'));
normA = 1;
normb = 1;
normR = 1;

C = Corig / normC;

scale = struct();
scale.normA = normA; scale.normb = normb;
scale.normR = normR; scale.normC = normC;
optSDPDAL.scale = scale;

% tolerance for SDPDAL
tol = tol_max;

% string format (for logging)
str_header = 'iter       obj          cut        nrmG    rank    lambda   ratio-1  time\n';
str_header_l = '-----------------------------------------------------------------------------------\n';
str_data = '%4d %12.6e % 12.6e %6.2e %5d   %6.2e %6.2e  %.2f';

tt = tic;
% main loop
for iter = 1:opts.maxit
    % prepare f
    f = mc_factory(opts.entropy, C, alpha, lambda);
    
    % call ALM_SSN
    optSDPDAL.R0 = x;
    optSDPDAL.tol = tol;
    [x,~,~,~,out] = ALM_SSN(n,[],[],AIop,bI,'oblique',f,[],p,optSDPDAL);
    p = size(x, 1);

    % now check the rank of x
    [~, D, V] = svd(x, 'econ');
    d = diag(D);
    rank_x = sum(d / d(1) > 1e-5);
    new_p = max(3, min(p, ceil(1.1 * rank_x)));
    if new_p < p
        x = D(1:new_p,1:new_p) * V(:,1:new_p)';
    end
    
    % compute ratio
    ratio = sum(d) / d(1);
    
    if opts.verbose
        if iter == 1
            fprintf(str_header);
            fprintf(str_header_l);
        end
        %cut = fun(x, Corig);
        [cut, ~] = rounding_sdp_gw(Corig, x, 50);
        fprintf(str_data, iter, out.pobj, cut, out.snrmG, p, lambda, ...
            ratio - 1, toc(tt));
        if new_p < p
            fprintf('  p:[%d -> %d]', p, new_p);
        end
        fprintf('\n');
    end
    
    % termination rule
    if ratio - 1 < opts.tol || strcmp(opts.entropy, 'none')
        cut = round(out.pobj);
        break;
    end
    
    % update sub_solver
    if ratio - 1 < 1e3 * opts.tol
        optSDPDAL.sub_solver = 1;
    end
    
    % update tol
    tol = max(tol_min, tol * 0.5);
    
    % update lambda and p
    lambda = max(1.1, min(5, ratio ^ rho)) * lambda;
    if new_p < p
        p = new_p;
    end
end

ret = struct();
ret.time = toc(tt);
ret.R = x;
ret.iter = iter;
ret.obj = cut;
ret.pinf = norm(dot(x,x,1) - 1, inf);
ret.p = p;

end

function options = default_SDPDAL_options(opts)

% options
options = [];
options.verbosity = opts.verbose > 1;
options.sub_solver = 3;
options.max_iter = 500;
options.maxitersub = 500;
options.tau = 0.8;
options.rho = 1.1;
options.sigma0  = 0.1;
options.sigma_min = 1e-2;
options.ALM_step = 1;
options.debug = 1;
options.comp_eigS = 0;
options.use_lobpcg = true;
options.use_sdpnal = 0;
options.gtol_ratio0 = 1e-3;
options.flag_decrease_p = true;
options.escape_saddle = true;
options.crit = 3;
options.tol = 1e-6;

end

function [f,g] = fun(X,C)
    % max Tr(C*Z), s.t., Z_ii = 1, Z psd
    % low rank model:
    % Z = X'*X, X = [X_1, ..., X_n], X is a p by n matrix
    % max Tr(C*X'*X), s.t., ||X_i|| = 1,
    g = 2*(X*C);
    f = sum(dot(g,X))/2;
end
