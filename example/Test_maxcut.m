function Test_maxcut(with_cutting_planes)
% maxcut problem (with cutting plane)
% set with_cuttin_planes = true to test max-cut with cutting planes

if nargin < 1
    with_cutting_planes = false;
end

% file path
file_dir = '..';
addpath(genpath(strcat(file_dir,'/src')));
addpath(strcat(file_dir,'/ARNT/solvers'));
addpath(strcat(file_dir,'/ARNT/manifolds-manopt'));
data_dir = [file_dir,'/sdp_data/'];

% gset file list
src = [data_dir 'Gset']; % for G set
probnames = {'g01', 'g02', 'g03', 'g04', 'g05'};

table_str = '';

% load main data

for ind_names = 1:numel(probnames)
    name = probnames{ind_names};
    rank_thres = 1e-5;
    load([src '/' name], 'A', 'n');
    d = sum(A, 2);
    Corig = (A - spdiags(d, 0, n, n)) / 4;

    % load cut data
    nc = ceil(sqrt(n/2));
    load([src '/cuts/' name '_cuts_T' num2str(nc) '.mat'], 'Bt');

    if with_cutting_planes
        % adding cutting planes using Bt
        AIop = mc_gen_AIop(Bt);
        bI = ones(numel(Bt), 1); 
    else
        AIop = []; bI = [];
    end

    p = min(ceil(sqrt(2 * n)), 200);
    R0 = randn(p, n);
    R0 = R0 ./ sqrt(sum(R0 .^ 2, 1));

    % scale parameter
    % Ai has norm == sqrt(1.5)
    % diag(X) = e
    % C has norm <= 1
    %normC = 1;
    normC = max(1, norm(Corig, 'fro'));
    normA = 1;
    normb = 1;
    normR = 1;

    C = Corig / normC;

    scale.normA = normA; scale.normb = normb;
    scale.normR = normR; scale.normC = normC;

    % construct Aop
    Aop = [];
    b = [];

    % construct f
    f = struct();
    f.obj_grad = @mc_obj_grad;
    f.hess = @(R, U, data, C) 0;
    f.data = {C};

    % options
    opts = [];
    opts.verbosity = 1;
    %opts.record_file = ['1.log'];
    %opts.record_file = [basename '.log'];
    opts.sub_solver = 3; % 1: arnt  3: RGBB  0:adpative

    opts.R0 = R0;
    opts.max_iter = 500;
    opts.tau = 0.8;
    opts.rho = 1.1;
    opts.sigma0  = 0.1;
    opts.sigma_min = 1e-2;
    opts.nu_min = 1e-2;
    opts.ALM_step = 1;
    opts.debug = 1;
    opts.comp_eigS = 1;
    opts.use_lobpcg = true;
    opts.adjust_p = 1;
    opts.use_sdpnal = 0;
    opts.gtol_ratio0 = 1e-3;
    %opts.gtol_ratio0 = 2;
    opts.flag_decrease_p = true;
    opts.favor_deltak = true;
    opts.disable_auto_sigma0 = 1;
    %opts.scale_C = scale_C;
    opts.scale = scale;
    opts.tol = 1e-6;
    opts.escape_saddle = true;
    opts.early_stop = false;
    opts.crit = 0;

    % call SDPDAL
    [R_,y,S_,Z_,out]=SDPDAL(n,Aop,b,AIop,bI,'oblique',f,[],p,opts);
    % rounding
    [cut, ~] = rounding_sdp_gw(-Corig, R_, 50);

    X = {};
    X{1} = (R_'*R_);
    S = {};
    S{1} = S_;
    Z = {};
    Z{1} = Z_;

    iter = out.iter;
    time = out.time;
    pinf = out.deltak;
    dinf = 0;
    pobj = out.pobj;
    dobj = out.dobj;
    relgap = out.relgap;

    X_item = X{1};
    [n,~] = size(X_item);
    m = length(y);



    X_eig = eig(X_item);
    S_eig = eig(S{1});


    X_rank = sub_comp_rank(X_eig,struct('thres',rank_thres));
    S_rank = sub_comp_rank(S_eig,struct('thres',rank_thres));

    etaK2 = out.etaK2;
    etaC1 = out.deltaC; avginit = out.HessianVector;
    sub_iter = out.sub_iter;
    etaK1 = 0;
    etaC2 = out.deltaX;
    etaC3 = out.deltaC2;

    p = out.p;
    marker = '';
    if out.flag >= 2; marker = '*'; end

    table_str = [table_str marker name];
    table_str = [table_str, sprintf(' & %d & %14.5e & %6d & %8.1e & %8.1e & %8.1e & %8.1e & %8.1e & %4d & %.1f & %d ', ...
        n, pobj, cut, pinf, relgap, etaK2, etaC1, etaC3, iter, time, X_rank)];
    table_str = [table_str ' \\ \hline' newline];

end

disp(newline);
disp(table_str);



end
