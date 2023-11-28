function [x, f, out] = PGBB(x, fun, proj, opts, varargin)

% Riemannian gradient method with BB step size
%   min F(x), s.t., x in M
%
% Input:
%           x --- initial guess
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [X, out]= OptStiefelGBB(X0, @fun, opts, data1, data2);
%           M --- projection operator
%
%        opts --- option structure with fields:
%                 record = 0, no print out
%                 maxit       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
%                 alpha       initial step size
%        rhols, eta, nt       parameters in line search
%   
% Output:
%           x --- solution
%           f --- function value at x
%         out --- output information
% -----------------------------------------------------------------------
% Reference: 
%  J. Hu, A. Milzark, Z. Wen and Y. Yuan
%  Adaptive Quadratically Regularized Newton Method for Riemannian Optimization
%
% Author: J. Hu, Z. Wen
%  Version 1.0 .... 2017/8

if nargin < 3
    error('at least three inputs: [x, f, out] = RGBB(x, fun, M, opts)');
elseif nargin < 4
    opts = [];
end

% termination rule
if ~isfield(opts, 'gtol');      opts.gtol = 1e-6;  end % 1e-5
if ~isfield(opts, 'xtol');      opts.xtol = 1e-6;  end % 1e-6
if ~isfield(opts, 'ftol');      opts.ftol = 1e-13; end % 1e-13

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'alpha');     opts.alpha  = 1e-3;   end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-6;   end
if ~isfield(opts, 'eta');       opts.eta  = 0.2;      end
if ~isfield(opts, 'gamma');     opts.gamma  = 0.85;   end
if ~isfield(opts, 'STPEPS');    opts.STPEPS  = 1e-10; end
if ~isfield(opts, 'nt');        opts.nt  = 3;         end % 3
if ~isfield(opts, 'maxit');     opts.maxit  = 200;   end
if ~isfield(opts, 'eps');       opts.eps = 1e-14;     end
if ~isfield(opts, 'record');    opts.record = 0;      end
if ~isfield(opts, 'radius');    opts.radius = 1;      end
if isfield(opts,  'nt');         opts.nt = 5;          end

fid = 1;
if opts.record
    if isfield(opts, 'record_fid')
        fid = opts.record_fid;
    elseif isfield(opts, 'recordFile')
        fid = fopen(opts.recordFile,'w+');
    end
end

% initial guess


% copy parameters
gtol = opts.gtol;
xtol = opts.xtol;
ftol = opts.ftol;
maxit = opts.maxit;
rhols = opts.rhols;
eta   = opts.eta;
gamma = opts.gamma;
record = opts.record;
nt = opts.nt;
alpha = opts.alpha;

% initial function value and gradient
[f,g] = feval(fun, x, varargin{:});

d = proj(x - g) - x;

nrmG = norm(g,'fro');




% initial iter. information 
out.nfe = 1; Q = 1; Cval = f; out.fval0 = f;

%% Print iteration header if debug == 1
if opts.record
    fprintf(fid,'\n%6s %15s %15s  %16s %9s %9s %5s %6s\n', ...
        'Iter', 'f(X)', 'Cval', 'nrmG', 'XDiff', 'FDiff', 'nls', 'alpha');
end

if record == 10; out.fvec = f; end
out.msg = 'max_it'; 
data = '';

% loop
for iter = 1:maxit
    
    xp = x; gp = g; fp = f; dp = d;
   
    nls = 1; deriv = rhols*nrmG^2; 
    
    
    

    % curvilinear search
    while 1
        
        x = proj(xp + alpha * d);
        [f,g] = feval(fun, x, varargin{:});
        
        out.nfe = out.nfe + 1;
        if f <=  Cval - alpha*deriv || nls >= 3
            break
        end
        alpha = eta*alpha;
        nls = nls+1;
    end
    
    stop = proj(x - g) - x;
    nrmG = norm(stop,'fro')/(1+norm(x,'fro'));

    out.nrmGvec(iter) = nrmG;

    if record == 10; out.fvec = [out.fvec; f]; end
    
    % difference of x

    s = x - xp;
  
    
    XDiff = norm(s,'inf')/alpha; % (relative Xdiff) ~ g
    FDiff = abs(f-fp)/(abs(fp)+1);
    
    % ---- record ----
    if opts.record
        fprintf(fid,...
            '%6d %20.13e %20.13e %9.2e %9.2e %9.2e %2d %9.2e\n', ...
            iter, f, Cval, nrmG, XDiff, FDiff, nls, alpha);
    end
    
    
    % check stopping
    crit(iter) = FDiff;
    mcrit = mean(crit(iter-min(nt,iter)+1:iter));
    
    % ---- termination ----
    if nrmG < gtol %|| XDiff < xtol || FDiff < ftol
    %if iter > 100
        %     if nrmG < gtol || XDiff < xtol || mcrit < ftol
        %    if nrmG < gtol
        out.msg = 'cvg';
        if nrmG  < gtol, out.msg = strcat(out.msg,'_g'); end
        if XDiff < xtol, out.msg = strcat(out.msg,'_x'); end
        %         if FDiff < ftol, out.msg = strcat(out.msg,'_f'); end
        if mcrit < ftol, out.msg = strcat(out.msg,'_mf'); end
        break;
    end
    
    % difference of gradient
    d = proj(x - g) - x;
    y = d - dp;
   
    
    % BB step size
    sy = abs(iprod(s,y));
    if sy > 0
        if mod(iter,2)==0; alpha = norm(s, 'fro')^2/sy;
        else alpha = sy/ norm(y, 'fro')^2; end
        % safeguarding on alpha
        alpha = max(min(alpha, 1e20), 1e-20);
    end
    
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + f)/Q;
    %Cval = f;       
end

if opts.record && isfield(opts, 'recordFile')
    fclose(fid);
end
out.XDiff = XDiff;
out.FDiff = FDiff;
out.mcrit = mcrit;
out.nrmG = nrmG;
out.fval = f;
out.iter = iter;
out.avginit = 0;


% Euclidean inner product
function a = iprod(x,y)
  a = real(sum(sum(conj(x).*y)));
end

end
