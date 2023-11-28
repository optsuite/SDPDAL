function [S, Out] = findS(S0, U, opts)

n = numel(S0);
%--------------------------------------------------------------------------
if ~isfield(opts, 'xtol');      opts.xtol = 1e-12; end
if ~isfield(opts, 'gtol');      opts.gtol = 1e-3; end
if ~isfield(opts, 'ftol');      opts.ftol = 1e-14; end
if ~isfield(opts, 'rho1');      opts.rho1  = 1e-4; end
if ~isfield(opts, 'rho2');      opts.rho2  = 0.9; end
if ~isfield(opts, 'm');         opts.m  = 5; end
if ~isfield(opts, 'maxit');     opts.maxit  = 1000; end
if ~isfield(opts, 'eps');       opts.eps = 1e-6; end
if ~isfield(opts, 'record');    opts.record = 0; end
if ~isfield(opts,'itPrint');    opts.itPrint = 1;   end
if ~isfield(opts, 'nls');       opts.nls = 5; end

% copy parameters
xtol = opts.xtol;
ftol = opts.ftol;
gtol = opts.gtol;
maxit = opts.maxit;
m    = opts.m;
%--------------------------------------------------------------------------

% Initial function value and gradient
% prepare for iterations\
S = S0;
[f,  g] = FG(U, S0);
Out.f = [];  Out.nfe = 1;

nrmx = norm(S, 'fro');
nrmG0 = norm(g, 'fro');

% set up storage for L-BFGS
SK = zeros(n,m);		  % S stores the last ml changes in x
YK = zeros(n,m);		  % Y stores the last ml changes in gradient of f.
istore = 0; pos = 0;  status = 0;  perm = [];
rho = [];

% main loop
for iter = 1:maxit
    % begin line search
    % clear workls;
    % store old point
    Sp = S;   nrmxp = nrmx;
    fp = f;   gp = g;   %nrmGp = nrmG;

    % compute search direction
    % if the first iteration
    if istore == 0
        d = -g;
    else
        d = LBFGS_Hg_Loop(-g);
        if sum(isnan(d) > 0); d = -g; end
    end

    deriv = sum(sum(d .* g));
    normd = norm(d, 'fro');
    if deriv > 0
        d = -g; deriv = -g' * g;
    end
    
    stp = 1;
    S = Sp + d;
    [f, g] = FG(U, S);
    
    k_inner = 0;
    while(k_inner <= opts.nls && f > fp + opts.rho1*stp*deriv + 1.0e-6)
        k_inner = k_inner+1;
        stp = 0.2^k_inner;
        S = Sp + stp*d; % backtracking
        [f, g] = FG(U, S);
    end % loop for while
    Out.nfe = Out.nfe + k_inner+1;
  

    % s = x - xp = stp*d;  ==> ||s|| = stp*||d||
    nrms = stp*normd;
    % compute stopping
    diffX = nrms/max(nrmxp,1);

    % now, update normG
    nrmG =  norm(g, 'fro');
    Out.nrmG =  nrmG;
    Out.f = [Out.f; f];

    cstop = nrmG / nrmG0 < gtol || (diffX < xtol) && (abs(fp-f)/(abs(fp)+1)) < ftol;
    
    if cstop; break; end

    nrmx = norm(S);

    %----------------------------------------------------------------------
    % save for L-BFGS
    ygk = g-gp;		s = S-Sp;

    %Check to save s and y for L-BFGS.
    if sum(sum(ygk .* ygk)) > 1e-20
        istore = istore + 1;
        pos = mod(istore, m);
        if pos == 0; pos = m; end
        YK(:,pos) = reshape(ygk, n, 1);  SK(:,pos) = reshape(s, n , 1);
        rho(pos) = 1/(sum(sum(ygk .* s)));
        
        if istore <= m
            status = istore;
            perm = [perm, pos]; 
        else
            status = m; 
            perm = [perm(2:m), perm(1)];
        end
    end

end


Out.msg = 'converge';
Out.iter = iter;
Out.nge = Out.nfe;
Out.fval = f;


% computer y = H*v where H is L-BFGS matrix
    function y = LBFGS_Hg_Loop(dv)
        q = reshape(dv, n, 1);   alpha = zeros(status,1);
        for di = status:-1:1
            k = perm(di);
            alpha(di) = (q'*SK(:,k)) * rho(k);
            q = q - alpha(di)*YK(:,k);
        end
        y = q/(rho(pos)* (sum(sum(ygk .* ygk))));
        for di = 1:status
            k = perm(di);
            beta = rho(k)* (y'* YK(:,k));
            y = y + SK(:,k)*(alpha(di)-beta);
        end
        y = reshape(y, size(S0));
    end

end


function [f, g] = FG(U, S)
    b = ones(size(U, 1), 1);
    Z = U * S;
    diagZ = sum(Z .^ 2, 2);
    f = norm(diagZ - b) ^ 2 / 4;
    g = U' * (diag(diagZ - b) * Z);
end