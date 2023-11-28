clear;
rng(31415);
file_dir = '.';
root_dir = './..';
addpath(file_dir);
data_dir = [root_dir '/sdp_data/'];
addpath(strcat(file_dir,'/ARNT/solvers'));
addpath(strcat(file_dir,'/ARNT/manifolds-manopt')); 
addpath(genpath(strcat(file_dir,'/manopt/')));
addpath(genpath('mexfun'));

% nc_type
% 0: nc = n
% 1: nc = ceil(sqrt(2n))
% 2: nc = ceil(sqrt(n/2))
nc_type = 2;

% ep_type|alpha|lambda
% 'N': no entropy
% 'T'|'R'
ep_type = 'T';
alpha = 2;
lambda = 10;

src = '../sdp_data/Gset'; % for G set
file_list = dir([src '/g*.mat']);
nfile = numel(file_list);
[~, ~] = mkdir(src, 'cuts');
for i = 1:nfile
    [~, name, ~] = fileparts(file_list(i).name);

    load([src '/' name], 'A', 'n');
    d = sum(A, 2);
    C = (A - spdiags(d, 0, n, n)) / 4;
    p = min(ceil(sqrt(2 * n) / 2), 200);
    R0 = randn(p, n);
    R0 = R0 ./ sqrt(sum(R0 .^ 2, 1));
    
    % construct Aop
    Aop = [];
    b = [];
    
    % construct f
    f = mc_factory(ep_type, C, alpha, lambda);
    
    
    opts = [];
    opts.verbosity = 2;
    %opts.record_file = [basename '.log'];
    opts.sigma_min = 0.1;
    opts.sub_solver = 1;
    opts.max_iter = 100;
    opts.tau = 0.8;
    opts.rho = 2;
    opts.sigma0 = 1e4;
    opts.nu0 = 1e4;
    opts.use_mex = 1;
    opts.ALM_step = 1.618;
    opts.debug = 1;
    opts.comp_eigS = 1;
    opts.adjust_p = 1;
    opts.use_sdpnal = 0;
    opts.gtol_ratio0 = 1e2;
    opts.flag_decrease_p = true;
    opts.favor_deltak = true;
    opts.R0 = R0;
    opts.crit = 3;
    %opts.disable_auto_sigma0 = true;
    
    %[R_,y_,S_,ret]=ALM_SSN(A_sub,b_sub,C_sub,p,opts);
    [R_,y_,S_,ret]=ALM_SSN(n,Aop,b,[],[],'oblique',f,[],p,opts);
    X = R_' * R_;
    
    fprintf('generating cuts ...\n');
    if nc_type == 0
        nc = n;
    elseif nc_type == 1
        nc = ceil(sqrt(2 * n));
    elseif nc_type == 2
        nc = ceil(sqrt(n / 2));
    else
        error('unknown nc_type.');
    end
    [pos, sgn] = mex_cutting_plane_mc(X, nc);
    
    fprintf('generating Bt ...\n');
    Bt = cell(1, nc);
    for ii = 1:nc
        data = [pos(1, ii), pos(2, ii), sgn(1, ii); ...
                pos(1, ii), pos(3, ii), sgn(2, ii); ...
                pos(2, ii), pos(3, ii), sgn(3, ii); ...
                n, n, 0];
        Bt{ii} = spconvert(data);
    end
    outfile = [src '/cuts/' name '_cuts_' ep_type num2str(nc) '.mat'];
    save(outfile, 'Bt', '-v7.3');
end