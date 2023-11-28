function [R,y,S,Z,ret]=SDPDAL(n, A, b, AI, bI, manifold, f, h, p, opts)
%-----------------------------------------------------------------------
% Consider the optimization problems on manifold:   
%    $$\min  f(X) + h(X)  s.t. X \in D, A(X) = b, AI(X) <= bI, $$   
% where $D = \{X: B(X)=b_0, X\succeq 0\}$, which has the manifold 
% structure by factorizing $X = R^TR$. 
%
% Input:
%           n --- dimension 
%           A --- linear equality constraints with fileds
%                  A.applyA(R)      A(R^TR)
%                  A.applyAT(y)     A^T(y)
%                  A.applyAUR(U,R)  A(R^TU + U^TR) 
%           b --- right hand side of linear equality constraints  
%          AI --- linear inequality constraints  
%          bI --- right hand side of linear inequality constraints
%    manifold --- manifold structures 
%          f  --- function value, gradient and Hessian  of f(R^TR) with respect to R                
%          h  --- function value, proximal operator and its Jacobian of nonsmooth term
%          p  --- estimation of rank
%
%        opts --- options structure with fields
%                  max_iter    max number of outer iterations
%                  maxitersub  max number of inner iterations
%                  sigma0      initial guess of penalty paramter for equality and inequality constraints 
%                  nu0         initial guess of penalty paramter for nonsmooth term h
%                  tol         accuracy tolerance for solving the SDP problem.
%                  sub_solver  subproblem solver  0: adaptive selection 1: ARNT 2: manopt trustregion 3: RGBB
% Output:
%           R --- primal variable 
%           y --- dual variable
%           S --- dual variable
%           Z --- dual variable 
%         ret --- output information  
%  
% Author: Haoyang Liu, Kangkang Deng, Zaiwen Wen
%   Version 0.1 .... 2023/06
%-------------------------------------------------------------------------
if nargin<8; h=[]; end
if nargin<9; p=0; end
if nargin<10; opts = []; end
if ~isfield(opts, 'max_iter'); opts.max_iter = 1e3; end
if ~isfield(opts, 'maxitersub'); opts.maxitersub = 100; end
if ~isfield(opts, 'rho'); opts.rho = 2; end
if ~isfield(opts, 'sigma0'); opts.sigma0 = 10; end
if ~isfield(opts, 'nu0'); opts.nu0 = 10; end
if ~isfield(opts, 'sigma_max'); opts.sigma_max = opts.sigma0*1e3; end
if ~isfield(opts, 'sigma_min'); opts.sigma_min = 1e-2; end
if ~isfield(opts, 'nu_max'); opts.nu_max = opts.nu0*1e3; end
if ~isfield(opts, 'nu_min'); opts.nu_min = 1e-2; end 
if ~isfield(opts, 'tol'); opts.tol = 1e-6; end
if ~isfield(opts, 'gtol0'); opts.gtol0 = 1e-2; end
if ~isfield(opts, 'tau'); opts.tau = 0.9; end
if ~isfield(opts, 'verbosity'); opts.verbosity = 0; end
if ~isfield(opts, 'ALM_step'); opts.ALM_step = 1.618; end
if ~isfield(opts, 'sub_solver');opts.sub_solver = 0; end
% 0: adaptive selection 1: ARNT 2: manopt trustregion 3: RGBB
if ~isfield(opts, 'debug'); opts.debug = 0; end
if ~isfield(opts, 'rankR_use_gap'); opts.rankR_use_gap = true; end
if ~isfield(opts, 'flag_increase_p'); opts.flag_increase_p = false; end
if ~isfield(opts, 'flag_decrease_p'); opts.flag_decrease_p = false; end
if ~isfield(opts, 'comp_eigS'); opts.comp_eigS = 1; end
if ~isfield(opts, 'pd_ratio'); opts.pd_ratio = 5; end
if ~isfield(opts, 'subtol_scheme'); opts.subtol_scheme = 0; end
% 0: gtol_ratio * sqrt(errA); 1: gtol_ratio * errA // adaptive gtol_ratio
if ~isfield(opts, 'gtol_ratio0'); opts.gtol_ratio0 = 1e3; end
if ~isfield(opts, 'record_file'); opts.record_file = ''; end
if ~isfield(opts, 'favor_deltak'); opts.favor_deltak = false; end
if ~isfield(opts, 'disable_cache'); opts.disable_cache = false; end
if ~isfield(opts, 'scale_C'); opts.scale_C = 1; end
if ~isfield(opts, 'escape_saddle'); opts.escape_saddle = true; end
if ~isfield(opts, 'escape_saddle_step'); opts.escape_saddle_step = 0; end
if ~isfield(opts, 'increase_p'); opts.increase_p = false; end
if ~isfield(opts, 'decrease_p'); opts.decrease_p = false; end
if ~isfield(opts, 'early_stop'); opts.early_stop = true; end
if ~isfield(opts, 'scale'); opts.scale = []; end
if ~isfield(opts, 'crit'); opts.crit = 0; end
% 0: standard; 0.5: standard (relax relgap)
% 1: early stop;  2: early stop + ignore deltaZ
% 3: early stop + ignore relgap; 4: early stop + ignore relgap + etaK2
if ~isfield(opts, 'par_sw'); opts.par_sw = default_par_sw; end
if ~isfield(opts, 'use_lobpcg'); opts.use_lobpcg = false; end
if ~isfield(opts, 'gtol_decrease'); opts.gtol_decrease = 0.8; end
if ~isfield(opts, 'gtol_min'); opts.gtol_min = 1e-8; end


prefix = 'log';
tol = opts.tol;
rho = opts.rho;
scale_C = opts.scale_C;
sub_solver = opts.sub_solver;
if sub_solver == 0
    sub_solver = 3;
end

%construct scale parameter
scale = opts.scale;
if isempty(scale)
    normA =  1; normb = 1; normR = 1;normC = 1;
else
    normA = scale.normA;
    normb = scale.normb;
    normR = scale.normR;
    normC = scale.normC;
end
    
    

% construct Aop
if isempty(A) % A is empty, which means no general linear equality constraints
    Aop = struct();
    Aop.applyA = @(varargin) 0;
    Aop.applyAUR  = @(varargin) 0;
    zero_h = @(V) 0;
    Aop.applyAT = @(y) zero_h;
    % set b = 0 locally
    b = 0;
    dm = 0;
elseif isnumeric(A)
    Aop = struct();
    Aop.applyA = @(varargin) applyA_numeric(A, varargin{:});
    Aop.applyAT = @(y) applyAT_numeric(A, y);
    dm = length(b);
elseif isstruct(A) % A is struct, only validate its usability
    Aop = A;
    if ~isfield(A, 'applyA')
        error('A.applyA is not defined.');
    end
    if ~isfield(A, 'applyAT')
        error('A.applyAT is not defined.');
    end
    dm = length(b);
else
    error('unsupported input of A.');
end

% construct AIop
if isempty(AI) % AI is empty, which means no general linear inequality constraints
    AIop = struct();
    AIop.applyA = @(varargin) 0;
    AIop.applyAUR  = @(varargin) 0;
    zero_h = @(V) 0;
    AIop.applyAT = @(y) zero_h;
    % set bI = 0 locally
    bI = 0;
    dmI = 0;
elseif isnumeric(AI)
    AIop = struct();
    AIop.applyA = @(varargin) applyA_numeric(AI, varargin{:});
    AIop.applyAT = @(y) applyAT_numeric(AI, y);
    dmI = length(bI);
elseif isstruct(AI) % AI is struct, only validate its usability
    AIop = AI;
    if ~isfield(AI, 'applyA')
        error('AI.applyA is not defined.');
    end
    if ~isfield(AI, 'applyAT')
        error('AI.applyAT is not defined.');
    end
    dmI = length(bI);
else
    error('unsupported input of AI.');
end

% check f (structure of function handle)
% [f, g_h] = f.obj_grad(R) returns the objective and the function handle of g
% f.hess(R, U) returns the p x n matrix R f''(R'R)[R'U + U'R]
if ~isfield(f, 'obj_grad')
    error('f.obj_grad is not defined.');
end
if ~isfield(f, 'hess')
    % currently we are using a second-order method to solve the
    % subproblem, thus f.hess must be valid
    % one may use first-order method in the future
    error('f.hess is not defined.');
end
if ~isfield(f, 'data')
    f.data = {};
end

% check manifold type
if strcmp(manifold, 'sphere')
    man_handle = @spherefactory;
elseif strcmp(manifold, 'oblique')
    man_handle = @obliquefactory;
elseif strcmp(manifold, 'stiefel')
    man_handle = @stiefelfactory_shape;
elseif strcmp(manifold, 'oblique_shape')
    man_handle = @obliquefactory_shape;
else
    error('unsupported manifold type %s.', manifold);
end

% check h (structure of function handle)
% [f, prox_h, prox_h_norm] = h.obj_prox(R, Z, nuk, data)
% h.hess(R, Z, nuk, data) returns R h''(R'R - Z/nuk)[U'R + R'U]
% different from f, h can be empty (Zero)
if isempty(h)
    h.obj_prox = @zero_obj_prox;
    h.hess = @(R, Z, nuk, U, data) 0;
    h.is_empty = true;
else
    if ~isfield(h, 'obj_prox')
        error('h.obj_prox is not defined.');
    end
    if ~isfield(h, 'hess')
        % currently we are using a second-order method to solve the
        % subproblem, thus h.hess must be valid
        % one may use first-order method in the future
        error('h.hess is not defined.');
    end
    h.is_empty = false;
end
if ~isfield(h, 'data')
    h.data = {};
end

if p == 0
    p = min(ceil(sqrt(2 * n)), 200);
end

R0 = randn(p,n); nrmR0 = dot(R0,R0,1);
R0 = bsxfun(@rdivide, R0, sqrt(nrmR0));
R0 = R0/norm(R0,'fro');
if isfield(opts, 'R0')
    R0 = opts.R0;
end
ftol = 1e-6;

optARNT.record = opts.verbosity > 1;
gtol_init = 1e-6;
optARNT.gtol = gtol_init;
optARNT.xtol = 1e-5;
optARNT.ftol = 1e-8;
M = man_handle(p,n,normR);
optARNT.solver_init = @RGBB;
optARNT.opts_init.record = opts.verbosity > 1;
optARNT.opts_init.tau   = 1e-6;%1e-3
optARNT.opts_init.maxit = 100;%2000;
optARNT.opts_init.gtol  = optARNT.gtol*1e3;
optARNT.opts_init.xtol  = optARNT.xtol*1e2;
optARNT.opts_init.ftol  = optARNT.ftol*1e2;
optARNT.solver_sub  = @RNewton;
optARNT.opts_sub.record = 0;
optARNT.opts_sub.tau    = 1e-6;
optARNT.opts_sub.maxitvec  = [100,110,120,150,150];
%optARNT.opts_sub.maxitvec  = [20,25,30,30,30];
optARNT.opts_sub.gtol   = optARNT.gtol*1e0;
optARNT.opts_sub.xtol   = optARNT.xtol*1e0;
optARNT.opts_sub.ftol   = optARNT.ftol*1e0;
optARNT.fun_TR = @fun_sub;
optARNT.fun_extra = @(data) data;
% dev-mc
%optARNT.tau = 1;
% dev-ncm
optARNT.tau = 0.01;
optARNT.usenumstab = 0;

% options for sub_comp_rank
opt_comp_rank = [];
if ~opts.rankR_use_gap
    opt_comp_rank = struct('gap', -1, 'thres', 1e-6);
end

sigmak = opts.sigma0;
nuk = opts.nu0;
maxitersub = opts.maxitersub;

cstop = 0;
t = tic;
m = length(b);
mI = length(bI);
y = zeros(m, 1);
yI = zeros(mI, 1);
if h.is_empty; Z = 0; else; Z = zeros(n); end
R = R0;

flag =0;
iter = 0;

eig_iter = 1e3;

% for lobpcg
if opts.use_lobpcg
    work_lpcg = struct();
    work_lpcg.ncv = p; % initial ncv
    work_lpcg.X0 = randn(n, p);
end

fid = 1;
if opts.verbosity > 0
    str0 = '    %6s';
    str1 = '      %6s';
    str2 = '        %6s';
    str3 = '    %12s';
    stra = ['\n%6s',str2,str3,str2, str0, str0, '    %3s', '  %6s','   %6s','      %6s','    %6s','   %6s', '    %6s', '     %6s','      %6s', str0, str0];
    str_head = sprintf(stra,...
        'iter','pobj',  'dobj', 'deltaC','relgap','etaK2','p', 'rank' ,'deltak','deltaZ','sigmak', 'nuk', 'siter/sinn','snrmG', 'etime', 'smsg');
    str_head_debug = sprintf('    %10s','  gtol_tgt');
    str1 = '  %3.2e';
    str_num = ['\n  %4d','  %+15.6e','  %+15.6e','  %+7.2e','  %7.2e','  %+7.2e','  %4d','  %4d','   %4.2e', '   %9.2e'...
        '  %4.2e', '  %4.2e', '   %4d',' %4d','      %8.2e','   %6.2f','    %-12s'];
    str_debug = ['   %4.2e'];
    
    if ~isempty(opts.record_file)
        if ~exist(prefix, 'dir')
            mkdir(prefix);
        end
        record_fname = [prefix '/' opts.record_file];
        fid = fopen(record_fname, 'w');
    end
end

p_flag = 0;

adj_sigma_flag = 0;
bad_etaK2 = 0;
relgap_cnt = 0;
deltak_cnt = 0;
arnt_use_rgbb_cnt = 0;
rgbb_sw_cnt = 0;
arnt_sw_cnt = 0;

%ftol_inc_step = 0;
rankR_same = 0;
rankR_reach_max = 0;
rankR_reach_max_tol = 2;
rankR = p;
rankR0_p = 0;

etaK2 = -1;
deltaC = -1;
errA = 1;
relgap = 1;

rankR_start = 0;
gtol_ratio = opts.gtol_ratio0;
deltak_reach_tol = false;
deltakZ_reach_tol = false;
RGBB_extra_iter = 0;

ARb = Aop.applyA(R) - b;
AIRb = AIop.applyA(R) - bI;
AIRb_plus = max(AIRb, 0);
deltak = normb * (norm(ARb) + norm(AIRb_plus))/(1+normb * (norm(b) + norm(bI)));

[~, ~, deltaZ] = h.obj_prox(R, 0, nuk, [], h.data{:});


out.nrmG = 1;
sub_iter = 0;
avginit = 0;
etaK2_cnt = 0;
etaK2_cnt_subopt = 0;
% dev-mc
% gtol_bnd = 0.1;
gtol_bnd = opts.gtol0;
escape_saddle_cnt = 0;
increase_p_cnt = 1;
etaK2_list = zeros(1, opts.max_iter);

% initialize output
ret = struct();
ret.flag = 99;
ret.msg = 'exceed max iteration';
ret.hist = struct();
ret.hist.deltak  = zeros(1, opts.max_iter);
ret.hist.deltaZ  = zeros(1, opts.max_iter);
ret.hist.deltaC  = zeros(1, opts.max_iter);
ret.hist.deltaC2 = zeros(1, opts.max_iter);
ret.hist.etaK2   = zeros(1, opts.max_iter);
ret.hist.relgap  = zeros(1, opts.max_iter);
ret.time_ARNT = 0;
ret.time_RGBB = 0;
ret.time_eigS = 0;

if dm == 0 && dmI == 0 && h.is_empty 
        gtol_decrease = 0.2;
else
    gtol_decrease = opts.gtol_decrease;
end
    
while iter<opts.max_iter && ~cstop
    iter = iter+1;
    t_sub = tic;
    R0 = R;
    ALM_step = opts.ALM_step;
    sub_solver_success = false;
    
    % etaK2 is bad for X consecutive iterations
    if bad_etaK2 > 0 && deltak_reach_tol
        % increase accuracy by X%
        gtol_ratio = max(gtol_ratio / 1.15, tol);
        % increase maxit for ARNT
        %maxitersub = maxitersub + 20;
        % temporarily decrease ALM_step
        if deltak_reach_tol
            ALM_step = 1;
        end
        % clear bad_etaK2
        bad_etaK2 = 0;
    end
    
    % set gtol
    gtol_min = 1e-8;
    if opts.subtol_scheme == 0
        gtol_ratio = opts.gtol_ratio0;
        if errA < tol && etaK2 < 5 * tol && opts.crit == 0
            gtol_ratio = opts.gtol_ratio0 / 5;
        end
        if errA < tol && etaK2 < 100 * tol && opts.crit == 1
            gtol_ratio = opts.gtol_ratio0 / 5;
        end
        gtol = max([gtol_ratio * sqrt(deltak), gtol_bnd, gtol_min*scale_C]);
    elseif opts.subtol_scheme == 1
        gtol_ratio = gtol_ratio / 1.02;
        gtol = max([min(1e-1, gtol_ratio * max(deltak, deltaZ)), gtol_min]);
        optARNT.opts_sub.maxitvec  = [40,40,40,40,40];
    else
        error('unknown subtol_scheme.\n');
    end

  
    
    
    % escape from saddle
    str_marker_s = '';
    if opts.escape_saddle && (etaK2 > 10 * errA || etaK2_stagnate_check(etaK2_list, iter, 3)) && errA < min(1e3 * tol, 1e-3) && rankR < p%size(R0, 1) && rankR/p<0.8
        if escape_saddle_cnt >= opts.escape_saddle_step %&& etaK2_stagnate_check(etaK2_list, iter, 4)
            [Q1, R1] = qr(R0);
            if ~exist('SS', 'var'); SS = []; end
            [RS, fl] = do_escape_saddle(S, Q1, R1, rankR, M, @fun_ARNT, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
            if fl == 0 % escape is success
                R0 = RS;
                str_marker_s = '[s]';
                escape_saddle_cnt = 0;
            else
                str_marker_s = '[/]';
            end
         else
             escape_saddle_cnt = escape_saddle_cnt + 1;
         end
    end
   
    
%     if(opts.decrease_p && iter>3 && etaK2< 10 * tol && deltak<5 * tol && rankR/p<0.7)
%        num_decrease = min(50, floor((p - rankR)/2)); 
%        [p,R0,M] = do_decrease_pp(R0,p,n,num_decrease,scale_C,man_handle);
%    end

    
    switch sub_solver
        case 1
            optARNT.gtol = gtol;%max(gtol_init,deltak);
            optARNT.ftol = 1e-12; %min(deltak*1e-3,ftol);
            optARNT.opts_init.gtol  = optARNT.gtol*1e2;% 1e1
            optARNT.opts_init.ftol  = optARNT.ftol*1e3;% 1e2
            optARNT.maxit = maxitersub;
            if ~isempty(opts.record_file)
                optARNT.record_fid = fid; %optARNT.record_file = 1;
            end
            optARNT.switch = true;
            optARNT.maxit_sw = 50 + RGBB_extra_iter;
            optARNT.do_init = true;
            optARNT.zeta = 0.2;
            
            if ~opts.disable_cache
                optARNT.fun_extra = @(data) prepare_cache_ARNT(data, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
            end
            
            % set grad and hessian
            optARNT.hess = @(R, U, data) hess_ARNT(R, U, data,Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
            optARNT.grad = @(R) grad_ARNT(R, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
            
           
            [R, ~, out] = arnt(R0, @fun_ARNT, M, optARNT, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
            if opts.sub_solver == 0
                if strcmp(out.msg(1:2), 'g:')
                    arnt_use_rgbb_cnt = arnt_use_rgbb_cnt + 1;
                else
                    arnt_use_rgbb_cnt = 0;
                end
                
                % arnt exit with 'g:' msg
                if arnt_use_rgbb_cnt > opts.par_sw.arnt_to_rgbb
                    arnt_use_rgbb_cnt = 0;
                    rgbb_sw_cnt = opts.par_sw.rgbb_sw_maxit; % do X gradient steps
                    sub_solver = 3;
                end
            end
            
            if out.nrmG > 2 * gtol
                RGBB_extra_iter = min(400, max(2 * RGBB_extra_iter, 10));
            else
                gtol_bnd = gtol_bnd * gtol_decrease;
            end
            sub_iter = sub_iter + out.iter;
            avginit = avginit + out.avginit;
            
            ret.time_ARNT = ret.time_ARNT + out.time_ARNT;
            ret.time_RGBB = ret.time_RGBB + out.time_RGBB + out.intime;
            
        case 2
            problem.M = M;
            data = struct();
            data.sigma = 0;
            
            %                 if use_mex
            %                     problem.cost = @(R) fun_only_with_mex(R, E,Ref_E,b,C,y,Z,lb,ub,sigmak,nuk);
            %                     problem.egrad = @(R) grad_with_mex(R, E,Ref_E,b,C,y,Z,lb,ub,sigmak,nuk);
            %                     problem.ehess = @(R, U) hess_with_mex(R, U, data,E, Ref_E,b,C,y,Z,lb,ub,sigmak,nuk);
            %                 else
            % 			        problem.cost = @(R) fun(R,A,b,C,y,sigmak);
            % 			        problem.egrad = @(R) grad(R,A,b,C,y,sigmak);
            % 			        problem.ehess = @(R, U) ehess(R, U,A,b,C,y,sigmak);
            %                 end
            options_manopt.verbosity = 0;
            options_manopt.maxiter = 100;
            options_manopt.tolgradnorm = gtol;
            [R, ~, out] = trustregions(problem, R0, options_manopt);
        case 3
           optRGB = optARNT;
            
           optRGB.xtol = 1e-5;  optRGB.ftol = 1e-8;  optRGB.gtol = gtol;
           optRGB.alpha = 1e-3; optRGB.rhols = 1e-6; optRGB.gamma = 0.85;
           optRGB.nt = 5;       optRGB.eta = 0.2;    optRGB.STPEPS = 1e-10;
           optRGB.maxit = maxitersub + RGBB_extra_iter;
           % dev-ncm
           %optRGB.maxit = 50;% + RGBB_extra_iter;
           optRGB.record = opts.verbosity > 1;
           optRGB.record_fid = fid;
           
           t_local = tic;
           [R, ~, out] = RGBB(R0, @fun_ARNT, M, optRGB, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
           t_local = toc(t_local);
           if opts.sub_solver == 0 && iter > 10
               if rgbb_sw_cnt > 0
                   rgbb_sw_cnt = rgbb_sw_cnt - 1;
               else
                   if out.iter ==  maxitersub + 400 && arnt_sw_cnt < opts.par_sw.sw_max
                       arnt_sw_cnt = arnt_sw_cnt + 1;
                       sub_solver = 1;
                   end
               end
           end
           if out.nrmG > gtol
               RGBB_extra_iter = min(400, max(2 * RGBB_extra_iter, 10));
           else
               gtol_bnd = gtol_bnd * gtol_decrease;
           end
           sub_iter = sub_iter + out.iter;
           avginit = avginit + out.avginit;
           
           ret.time_RGBB = ret.time_RGBB + t_local;
    end
    sub_time = toc(t_sub);
    acc_time = toc(t);
    
    % compute rank(R)
    rankR_p = rankR;
    [~,R1,~] = qr(R);
    diag_R1 = abs(diag(R1));
    rankR = sub_comp_rank(diag_R1, opt_comp_rank);
    
    RR = R' * R;
    ARb = Aop.applyA(R) - b;
    AIRb = AIop.applyA(R) - bI;
    AIRb_plus = max(AIRb, 0);
    
    deltak_p =deltak;
    deltak = normb * (norm(ARb) + norm(AIRb_plus))/(1+normb * (norm(b) + norm(bI)));
    %deltak = norm(ARb) / (1 + norm(b));
    
    deltak_ratio = deltak/deltak_p;
    
    y_p = y;
    sigmak_p = sigmak;
    nuk_p = nuk;
    
    % multiplier update: y, yI
    y = y-ALM_step*sigmak*ARb;
    yI = yI - ALM_step*sigmak*(AIRb + max(yI/sigmak - AIRb, 0));
    
    deltaZ_p = deltaZ;
    [~, W_h] = h.obj_prox(R, Z, nuk, [], h.data{:});
    dZ = W_h(speye(n))+ Z / nuk;
    deltaZ = norm(dZ, 'fro') / (max(1, norm(RR, 'fro'))+1);
    if h.is_empty
        deltaZ_ratio = 0;
    else
        deltaZ_ratio = deltaZ / deltaZ_p;
    end

    Z_p = Z;
%     plot(((eigs(Z,n))));
%     hold on
    % multiplier update: Z
    Z = Z-ALM_step*nuk*dZ;
    %Z = Z - 0.5*nuk*(RR-max(RR-Z/nuk,0));
    %Z =nuk*RR - max(RR,0);
    % if ~reach_deltak
    % 	sigtol = opts.tau;
    % end
    
    % if  (deltak_ratio < sigtol || deltak < opts.tol || sigmak == max_sigmak || adj_sigma_flag)
    % 	y = y-ALM_step*sigmak*ARb;
    % elseif iter<opts.max_iter-1
    % 	sigmak = min(opts.rho*sigmak,max_sigmak);
    % end
    
    % optimality condition: etaK2 (dual feas)
    
    % [~, g_h] = f.obj_grad(R, f.data{:});
    % [objR, gradR] = fun_ARNT(R,Aop, b, f, h, y, Z, sigmak, nuk);
    % objX = objR;
    % C = g_h(eye(n));
    
    % ATy_h = Aop.applyAT(y_p);
    % S1 = C - ATy_h(eye(n));
    % % FIX ME!! the formula of S is only correct on sphere manifold and oblique manifold
    % [u, BTu] = compute_u(manifold, gradR, R);
    % S = S1-BTu-Z;
    % % FIX ME!! this only works for the box constraint
    % dobj = sum(b.*y)+sum(u);
    
    % FIX ABOVE!!


    % apply AT (handle)
    ATy_h = Aop.applyAT(y);
    ATy = ATy_h(speye(n));
    AITy_h = AIop.applyAT(yI);
    AITy = AITy_h(speye(n));

    % call f
    [f1, g1_h] = f.obj_grad(R, f.data{:});

    % call h (compute RT(RRZ) )
    % f2 = h(prox(RRZ))
    % g2_h = handle of T(RRZ)
    [f2, ~, ~] = h.obj_prox(R, Z, nuk, [], h.data{:});

    objR = f1+f2;
    objX = objR;

    S1 = g1_h(speye(n))-ATy-AITy-Z;
    gradR = 2*R*S1;
    %gradR = grad_ARNT(R, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
    % FIX ME!! the formula of S is only correct on sphere manifold and oblique manifold
    [u, BTu] = compute_u(manifold, gradR, R, normR);
    S = S1-BTu;
    normS = norm(S, 'fro');

    if strcmp(manifold, 'sphere')
        dobj = sum(b.*y) + sum(bI.*yI) + normR^2 * u;
    elseif strcmp(manifold, 'oblique')
        dobj = sum(b.*y) + sum(bI.*yI) + normR^2 * sum(u);
    elseif strcmp(manifold, 'stiefel')
        dobj = sum(b.*y) + sum(bI.*yI) + normR^2 * trace(u);
    elseif strcmp(manifold, 'oblique_shape')
        dobj = sum(b.*y) + sum(bI.*yI) + normR^2 * sum(u);
    end
    if isfield(f,'conjugate') && isfield(h,'conjugate')  %% need to fix, if user not provide conjugate, we will not compute gap
        dobj = dobj - f.conjugate(ATy + BTu + Z + S, f.data{:}) -  h.conjugate(-Z,h.data{:});
    end
        
    % rescale     
    objR = objR*normC*normb/normA; 
    dobj = dobj*normC*normb/normA;
    % relative gap
    relgap = abs(objR - dobj)/(1+abs(objR)+abs(dobj));
    
    t_local = tic();
    if opts.comp_eigS > 0 && mod(iter - 1, opts.comp_eigS) == 0
        
        % compute etaK2
        etaK2_p = etaK2;
        if opts.use_lobpcg
            [V, S_eig] = lobpcg(work_lpcg.X0, S, tol);
            ncv = min(n, 2 + ceil(1.2 * sum(S_eig < 0)));
            if ncv > work_lpcg.ncv
                % padding
                work_lpcg.X0 = [V, randn(n, ncv - work_lpcg.ncv)];
            else
                work_lpcg.X0 = V(:, 1:ncv);
            end
            work_lpcg.ncv = ncv;
        else
            [S_eig] = eig(full(S), 'vector');
        end
        SS = [];
        %SS.V = V; SS.d = S_eig;
        
        S_eig_neg = S_eig(S_eig < 0);
        etaK2 = norm(S_eig_neg) / (1 + normS);
        %etaK2 = normC * norm(S_eig_neg) / (1 + normC * norm(S_eig));
        etaK2_list(iter) = etaK2;
        
        % adjust gtol
        etaK2_ratio = etaK2 / etaK2_p;
        if (deltak_ratio > opts.tau || deltak < tol) && etaK2_ratio > opts.tau && (deltak_reach_tol || opts.pd_ratio * deltak < etaK2)
       % if deltak<10*tol && 100 * deltak < etaK2 && etaK2_ratio > opts.tau
            % etaK2 is bad, increase bad_etaK2 count
            bad_etaK2 = bad_etaK2 + 1;
        else
            bad_etaK2 = 0;
        end
    end
    ret.time_eigS = ret.time_eigS + toc(t_local);
    
    % adjust sigmak such that deltak & etaK2 decrease at the same rate
    sigtol = opts.tau;
    if deltak_ratio > sigtol
        deltak_cnt = deltak_cnt + 1;
    else
        deltak_cnt = 0;
    end
%     if deltak > tol && deltak_cnt > 1 && adj_sigma_flag == 0 % deltak too large
%         sigmak = min(rho * sigmak, opts.sigma_max);
%         if opts.favor_deltak && (sigmak > sigmak_p && ~deltak_reach_tol || etaK2 < tol)
%             gtol_ratio = rho * gtol_ratio;
%         elseif deltak_reach_tol
%             gtol_ratio = gtol_ratio / rho;
%         end
%         adj_sigma_flag = 1;
%     elseif deltak < tol && etaK2 > tol*100 %&& opts.pd_ratio * deltak < etaK2 % etaK2 too large
%         sigmak = max(sigmak / sqrt(rho), opts.sigma_min);
%         adj_sigma_flag = 1;
%         deltak_reach_tol = true;
%     else
%         % relgap too large
%         if deltak < tol && etaK2 < tol && relgap > tol
%             relgap_cnt = relgap_cnt + 1;
%         else
%             relgap_cnt = 0;
%         end
%         if opts.favor_deltak || deltak_ratio > sigtol ^ 2
%             gtol_ratio = gtol_ratio / 1.02;
%         end
%         adj_sigma_flag = 0;
%     end
    
    if deltaZ_ratio > sigtol && deltaZ >= tol
        nuk = min(rho * nuk, opts.nu_max);
    elseif deltaZ < tol && etaK2 > tol*10
        nuk = max(nuk / rho, opts.nu_min);
    end

    if deltak_ratio > sigtol && deltak >= tol
        sigmak = min(rho * sigmak, opts.sigma_max);
    elseif deltak < tol && (etaK2 > tol*1.1 || deltaC > tol)
        sigmak = max(sigmak / rho, opts.sigma_min);
    end
    
    if max(deltak, deltaZ) < tol
        deltakZ_reach_tol = true;
    end
    
    if opts.subtol_scheme == 1 && (sigmak > sigmak_p || nuk > nuk_p) && ~deltakZ_reach_tol
        %gtol_ratio = gtol_ratio * rho;
    end
    
%     if sigmak > sigmak_p || nuk > nuk_p
%         gtol_ratio = gtol_ratio * sqrt(rho);
%     elseif sigmak < sigmak_p || nuk < nuk_p
%         gtol_ratio = gtol_ratio / sqrt(rho);
%     end
Use Control + Shift + m to toggle the tab key moving focus. Alternatively, use esc then tab to move to the next interactive element on the page.
Editing SDPDAL/src/SDPDAL.m at main · optsuite/SDPDAL

    
    
    %S = S1-BTu-Z;

    %RR = R'*R;
    scaleRR = RR*normb/normA;
    deltaC = normC*abs(sum(sum(scaleRR.*S)))/(1+norm(scaleRR,'fro') + normC*normS);
    deltaC2 = abs(sum(yI .* AIRb)) / (1 + norm(yI) + norm(AIRb));
    if opts.crit == 0.5 && max([deltak, deltaZ, deltaC, deltaC2, etaK2]) < tol && relgap < 2 * tol
        relgap_cnt = relgap_cnt + 1;
    end
    % set termination boolean
    %deltaX = norm(RR(RR<0),'fro')/(1+norm(RR,'fro'));
    cstop0 = (relgap < tol || relgap_cnt > 2 || opts.crit >= 3)...
        && deltak < tol && deltaC < tol ...
        && deltaC2 < tol && (deltaZ < tol || opts.crit == 2);
    cstop = cstop0 && (etaK2 < tol || opts.crit == 4);
    if cstop
        ret.flag = 0;
        ret.msg = 'converge';
    elseif cstop0 % etaK2 > tol
        if etaK2 < 5 * tol
            etaK2_cnt_subopt = etaK2_cnt_subopt + 1;
            etaK2_cnt = 0;
        else
            etaK2_cnt = etaK2_cnt + 1;
            etaK2_cnt_subopt = 0;
        end
        
        if etaK2_cnt_subopt > 2 && (opts.early_stop || opts.crit >= 1)
            cstop = true;
            ret.flag = 1; % sub optimal
            ret.msg = 'converge to sub-optimal';
        end
        
        if (opts.escape_saddle && etaK2_cnt > 100 || ~opts.escape_saddle && etaK2_cnt > 20) && (opts.early_stop || opts.crit >= 1)
            cstop = true;
            ret.flag = 2; % early stopping (saddle)
            ret.msg = 'potential saddle found';
        end
       %cstop = (relgap < tol || relgap_cnt > 2) && deltaZ<tol && etaK2 && deltak < tol && deltaC <tol && deltaX<tol;
    end
   % obj = opts.round(R)
    %if rankR == 1 || out.nrmG < tol  ||  (obj - opts.opt) < 1 || acc_time>opts.max_time cstop = 1; end
    % update errA (overall error w/o etaK2)
    errA = max([deltak, deltaC, deltaC2]);
    if opts.crit ~= 2; errA = max(errA, deltaZ); end
    if opts.crit ~= 3; errA = max(errA, relgap); end
    
    %cstop = out.nrmG<tol*1e1 && deltak<tol && deltaZ<tol;
    % print
    if opts.verbosity > 0
        % print header
        if iter == 1 || opts.verbosity > 1
            fprintf(fid, str_head);
            if opts.debug
                fprintf(fid, str_head_debug);
            end
        end
        
        % print iteration info
        switch sub_solver
            case 1
                fprintf(fid, str_num,iter,objR,dobj,max(deltaC,deltaC2), relgap, etaK2,p,rankR,deltak,deltaZ,sigmak_p,nuk_p, out.iter, out.avginit,out.nrmG, acc_time, out.msg);
                if opts.debug
                    fprintf(fid, str_debug, gtol);
                end
            case 2
                fprintf(fid, str_num,iter,objR,objX,dobj, relgap,etaK2,p,rankR, deltak,deltaZ,sigmak_p,out(end).iter, out(end).gradnorm, acc_time, '');
                if opts.debug
                    fprintf(fid, str_debug, gtol);
                end
            case 3
                fprintf(fid, str_num,iter,objR,dobj,max(deltaC,deltaC2), relgap, etaK2,p,rankR,deltak,deltaZ,sigmak_p,nuk_p, out.iter, out.avginit,out.nrmG, acc_time, out.msg);
                if opts.debug
                    fprintf(fid, str_debug, gtol);
                end
        end
        
    end
    
    % record iter history
    ret.hist.deltak(iter)  = deltak;
    ret.hist.deltaZ(iter)  = deltaZ;
    ret.hist.deltaC(iter)  = deltaC;
    ret.hist.deltaC2(iter) = deltaC2;
    ret.hist.etaK2(iter)   = etaK2;
    ret.hist.relgap(iter)  = relgap;
    
    % increase p
    if rankR==p
        rankR_reach_max = rankR_reach_max+1;
    else
        rankR_reach_max = 0;
    end
    
    if opts.flag_increase_p && rankR_reach_max >= rankR_reach_max_tol
        error('Not implemented'); % FIX ME!!
    end
    
    % decrease p
    if rankR - rankR_p <= 1
        if rankR_same == 0
            rankR_start = iter;
        end
        rankR_same = rankR_same+1;
    else
        rankR_same = 0;
    end
    
     if(opts.increase_p && iter>1 && etaK2> 5 * tol &&  p<n && rankR/p>0.8  && increase_p_cnt<8)% && mod(iter - 1, opts.comp_eigS) == 0
        increase_p_cnt = increase_p_cnt  + 1;
        S_eig = S_eig(S_eig<-1e-6); num_increase = length(S_eig);
        %[~, num_increase] = max(diff(S_eig)./abs(S_eig(2:end)));
        num_increase = min(num_increase,10);
        if num_increase > 1
        [p,R,M] = do_increase_pp(S,R,p,n,num_increase, scale_C,man_handle,@fun_ARNT,Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk);
        end
     end
    
    if opts.flag_decrease_p && ~cstop && rankR_same > 5 && (errA < 1e-3 || rankR_start > 5)
        if etaK2 > 1e-3
            p_extra = 15;
        elseif etaK2 > 1e-4
            p_extra = 10;
        elseif etaK2 > 1e-5
            p_extra = 7;
        else
            p_extra = 5;
        end
        new_p = min(n, max(p_extra + rankR, ceil(1.1 * rankR)));
        % only allow decreasing
        if new_p < p
            if opts.verbosity
                fprintf(fid, '[p: %d -> %d]', p, new_p);
            end
            p = new_p;
            [R, M] = do_decrease_p(new_p, R, normR, man_handle);
        end
    end
    
    if opts.verbosity
        fprintf(fid, '%s', str_marker_s);
    end
 
    
end

% toc here
tsolve_ALM = toc(t);

sub_iter = sub_iter/iter;

 gradR = 2*R*S1;

[u, BTu] = compute_u(manifold, gradR, R, normR);
u = u(:);
len_u = length(u);
y_original = zeros(dm+dmI+len_u, 1);
y_original(1:len_u) = u*normA;
if dm > 0
    y_original((len_u+1):(len_u+dm)) = y;
end
if dmI > 0
    y_original((len_u+dm+1):(len_u+dm+dmI)) = yI;
end

% re-scale solution varible
y = y_original*normC/normA;
R = sqrt(normb/normA)*R;
S = normC*S;
Z = normC*Z;

ret.time = tsolve_ALM;
ret.iter = iter;
ret.deltak = deltak;
ret.sigmak = sigmak;
ret.snrmG = out.nrmG;
ret.X = R'*R;
ret.Z = Z;
ret.BTu = BTu;
ret.u = u;
ret.gradR = gradR;
ret.W = R'*R-dZ;
ret.nu = nuk;
ret.pobj = objR;
ret.dobj = dobj;
ret.sub_iter = sub_iter;
ret.deltaX = sum(sum(scaleRR.*Z))/(1+norm(scaleRR,'fro') +norm(Z,'fro'));
ret.deltaC = deltaC;
ret.deltaC2 = deltaC2;
ret.p = p;
ret.rankR = rankR;
ret.HessianVector = avginit/iter;
ret.etaK2 = etaK2;
ret.relgap = relgap;
ret.hist.deltak  = ret.hist.deltak(1:iter);
ret.hist.deltaZ  = ret.hist.deltaZ(1:iter);
ret.hist.deltaC  = ret.hist.deltaC(1:iter);
ret.hist.deltaC2 = ret.hist.deltaC2(1:iter);
ret.hist.etaK2   = ret.hist.etaK2(1:iter);
ret.hist.relgap  = ret.hist.relgap(1:iter);
ret.nrmG = out.nrmG;

if opts.verbosity
    hrule = repmat('-', 1, 80);
    fprintf(fid, '\n%s\n', hrule);
    fprintf(fid, '- ALM_SSN OUTPUT\n');
    fprintf(fid, [hrule, '\n']);
    fprintf(fid, '  exit code = %d (%s)\n', ret.flag, ret.msg);
    fprintf(fid, [hrule, '\n']);
    fprintf(fid, '  iter = %d\n', iter);
    fprintf(fid, '  time = %.2f\n', ret.time);
    fprintf(fid, '    time (ARNT) = %.2f\n', ret.time_ARNT);
    fprintf(fid, '    time (RGBB) = %.2f\n', ret.time_RGBB);
    fprintf(fid, '    time (eigS) = %.2f\n', ret.time_eigS);
    fprintf(fid, '  pobj = %.12e,  dobj = %.12e, relgap = %.2e\n', ret.pobj, ret.dobj, relgap);
    fprintf(fid, '  feasibility:\n');
    fprintf(fid, '    deltak = %8.2e    deltaZ  = %8.2e\n', deltak, deltaZ);
    fprintf(fid, '    deltaC = %8.2e    deltaC2 = %8.2e\n', deltaC, deltaC2);
    fprintf(fid, '    etaK2  = %8.2e\n', etaK2);
    fprintf(fid, '  rank(R) = %d\n', rankR);
    fprintf(fid, [hrule, '\n']);
end
%fprintf(fid, '\nOUTPUT \n & %2d & %.4e & %.2e & %.2e & %.2f  & %.2e & %.2e & %.2e & %.2e & %2d \n',iter,objR,out.nrmG,relgap,tsolve_ALM, deltak, deltaZ, etaK2, deltaC, rankR );
if opts.verbosity > 0 && ~isempty(opts.record_file)
    fclose(fid);
end

end

function AX = applyA_numeric(A, R, U)

% for R'R
if nargin == 2
    RTR_vec = reshape(R'*R, [], 1);
    AX = A * RTR_vec;
end

% for U'R + R'U
if nargin == 3
    UR = U' * R;
    UR_vec = reshape(UR+UR', [], 1);
    AX = A * UR_vec;
end

end

% compute ATy in the form of handle ATy_h(R) = R * ATy
function ATy_h = applyAT_numeric(A, y)

n_sqr = size(A, 2);
n = round(sqrt(n_sqr));
ATy = reshape(A' * y, n, n);

% note: ATy is the form of handle
ATy_h = @(R) R * ATy;

end

function [f,g] = fun_ARNT(R, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk)

% apply A
ARb = Aop.applyA(R) - b - y/sigmak;
% apply AT (handle)
AARb_h = Aop.applyAT(ARb);

% apply AI
AIRb = AIop.applyA(R) - bI - yI/sigmak;
AIRb_plus = max(AIRb, 0);

% apply AIT (handle)
AAIRb_h = AIop.applyAT(AIRb_plus);

% call f
[f1, ~, g1] = f.obj_grad(R, f.data{:});

if isstruct(g1) && isfield(g1,'RR')
    data.RR = g1.RR; g1 = g1.grad;
else
    data = [];
end
    
% call h (compute RT(RRZ) )
% f2 = h(prox(RRZ))
% g2_h = handle of T(RRZ)
[f2, g2_h, g2_norm] = h.obj_prox(R, Z, nuk, data, h.data{:});

rg2h = 2*nuk*g2_h(R);
g = 2*g1+2*sigmak*(AARb_h(R) + AAIRb_h(R)) +  rg2h;

f = f1 + f2 + sigmak/2*(sum(AIRb_plus.^2) + sum(ARb.^2))+nuk/2*g2_norm^2;

end

% function f = fun_only_ARNT(R, Aop, b,C,y,Z,lb, ub, sigmak, nuk)
%
%     % apply A
%     ARb = Aop.applyA(R) - b - y/sigmak;
%
%     % box constraint
%     if ~isempty(lb) || ~isempty(ub)
%         RRZ = R' * R - Z / nuk;
%         HR = min(RRZ - lb, 0);
%     else
%         HR = 0;
%     end
%
%     g1 = 2*R*C;
%     f = sum(dot(g1,R))/2+sigmak/2*sum(ARb.^2)+nuk/2*norm(HR,'fro')^2;
%
% end

function [g,gg] = grad_ARNT(R, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk)

% apply A
ARb = Aop.applyA(R) - b - y/sigmak;
% apply AT (handle)
AARb_h = Aop.applyAT(ARb);

% apply AI
AIRb = AIop.applyA(R) - bI - yI/sigmak;
AIRb_plus = max(AIRb, 0);
% appli AIT (handle)
AAIRb_h = AIop.applyAT(AIRb_plus);

% call f
[~, ~, g1] = f.obj_grad(R, f.data{:});
if isstruct(g1) && isfield(g1,'RR')
    data.RR = g1.RR; g1 = g1.grad;
else
    data = [];
end
% call h (compute RT(RRZ) )
% f2 = h(prox(RRZ))
% g2_h = handle of T(RRZ)
[~, g2_h] = h.obj_prox(R, Z, nuk, data, h.data{:});

rg2h = 2*nuk*g2_h(R);
g = 2*g1+2*sigmak*(AARb_h(R) + AAIRb_h(R)) + rg2h;
if(nargout>1)
    gg = rg2h;
end
end

function hh = hess_ARNT(R, U, data, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk)

% read AARb from cache if possible
if isfield(data, 'AARb_h')
    AARb_h = data.AARb_h;
else % if not present, compute AARb_h
    % apply A
    ARb = Aop.applyA(R) - b - y/sigmak;
    % apply AT (handle)
    AARb_h = Aop.applyAT(ARb);
end

% read D/AAIRb from cache if possible
if isfield(data, 'D') && isfield(data, 'AAIRb_h')
    D = data.D;
    AAIRb_h = data.AAIRb_h;
else % if not present, compute D = AI(X) - bI - yI/sigma
    % apply AI
    AIRb = AIop.applyA(R) - bI - yI/sigmak;
    D = AIRb > 0;
    % apply AIT (handle)
    AAIRb_h = AIop.applyAT(max(AIRb, 0));
end

% apply A (on R'U + U'R)
AUR = Aop.applyAUR(R, U);
% apply AT (on AUR)
AAUR_h = Aop.applyAT(AUR);

% apply AI (on R'U + U'R)
AIUR = AIop.applyA(R, U);
% apply AIT (on D * AIUR)
AAIUR_h = AIop.applyAT(D .* AIUR);



% read g_h from cache if possible
if isfield(data, 'g_h')
    g_h = data.g_h;
else % if not present, compute g_h
    [~, g_h] = f.obj_grad(R, f.data{:});
end

% read prox_h from cache if possible
if isfield(data, 'prox_h')
    prox_h = data.prox_h;
else % if not present, compute prox_h
    [~, prox_h] = h.obj_prox(R, Z, nuk, data, h.data{:});
end

hh = 2*(g_h(U) +  sigmak*(AARb_h(U) + AAIRb_h(U))) + ...
    2*(f.hess(R, U, data, f.data{:}) + ...
         sigmak*(AAUR_h(R) + AAIUR_h(R))) + ...
         data.sigma*U + 2*nuk*prox_h(U) +2*nuk*h.hess(R, Z, nuk, U, data, h.data{:});




end

function data = prepare_cache_ARNT(data, Aop, b, AIop, bI, f, h, y, yI, Z, sigmak, nuk)

% read the current R
R = data.XP;

% call cache for f
if isfield(f, 'fun_extra')
    data = f.fun_extra(data);
end

if ~isfield(data, 'AARb_h')
    % apply A
    ARb = Aop.applyA(R) - b - y/sigmak;
    % apply AT (handle)
    AARb_h = Aop.applyAT(ARb);
    
    % save in cache
    data.AARb_h = AARb_h;
end

if ~isfield(data, 'AAIRb_h') || ~isfield(data, 'D')
    % apply AI
    AIRb = AIop.applyA(R) - bI - yI/sigmak;
    AIRb_plus = max(AIRb, 0);
    % apply AT (handle)
    AAIRb_h = AIop.applyAT(AIRb_plus);
    % compute D
    D = AIRb > 0;
    
    % save in cache
    data.AAIRb_h = AAIRb_h;
    data.D = D;
end

if ~isfield(data, 'g_h')
    % compute f.obj_grad
    [~, g_h] = f.obj_grad(R, f.data{:});
    
    % save
    data.g_h = g_h;
end

if ~isfield(data, 'RRZ') && ~h.is_empty
    data.RRZ = R' * R - Z / nuk;
end

if ~isfield(data, 'prox_h') && ~h.is_empty
    % compute h.obj_prox
    [~, prox_h] = h.obj_prox(R, Z, nuk, data, h.data{:});
    
    % save
    data.prox_h = prox_h;
end

end

function [f, prox_h, prox_h_norm] = zero_obj_prox(~, ~, ~, ~)
f = 0;
prox_h = @(V) 0;
prox_h_norm = 0;
end

function par = default_par_sw
par = struct();
par.arnt_to_rgbb = 2;
par.rgbb_sw_maxit = 8;
par.sw_max = 5;
end
