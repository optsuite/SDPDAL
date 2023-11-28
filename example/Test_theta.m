function Test_theta(test_type, test_id, rank_thres, compress)
%clear;

if nargin<1; test_type = 10; end
if nargin<2; test_id = 0;   end
if nargin<3; rank_thres = 1e-5;   end
if nargin<4; compress = true; end
% test_type
% 0: theta_new full
% 10: theta full
% 11: theta small
% 12: theta large
% 13: theta large (extra)
% -1: theta_new: brock200_1 -- brock200_4 (debug mode)
% -2: theta: 1xx.128
% -3: theta: 1dc.128, 1dc.256

% test_id
% 0: arnt only
% 1: grad
% 2: adaptive
% 10+x: rank-two formulation
% 99: SDPNAL+
if test_id >= 10 && test_id < 20
    variant = 1;
else
    variant = 0;
end
solver_id = mod(test_id, 10);

% file path
file_dir = '..';
root_dir = '../..';
addpath(file_dir);
addpath(genpath(strcat(file_dir,'/src')));
addpath(strcat(file_dir,'/ARNT/solvers'));
addpath(strcat(file_dir,'/ARNT/manifolds-manopt'));
data_dir = [file_dir,'/sdp_data/'];

    
% generate matnames


    src = [data_dir, 'theta']; % for theta
    save_root = strcat(file_dir,'/results/theta/');
    matnames = thetaprobs_modif;
    nfile = numel(matnames);



    prob_range = 1;


% preload rss_sol when variant == 1
if variant == 1
    load([src, '/mss_soln.mat'], 'mss_soln');
end

table_str = '';
comp_table_str = '';
detail_table_str = '';





for i = prob_range
    % load main data
    basename = matnames{i};
    fname = [src '/' basename '_Alt.mat'];
    prob_str = strrep(matnames{i}, '_', '-'); % for printing to latex
    
    
        % load graph
        graph = load(fname);
        E = sortrows(graph.E, 1);
        n = graph.n;
        
        num_E = size(E, 1);
        Ref_E = sparse(E(:, 1), E(:, 2), ones(num_E, 1), n,n)';
        
        % construct Aop
        Aop = struct();
        ratio = nnz(Ref_E) / (n * (n - 1) / 2);
        if ratio < 0.05 % sparse
            Aop.applyA = @(R) theta_applyA(E, R);
            Aop.applyAUR = @(R, U) theta_applyA(E, R, U);
            Aop.applyAT = @(y) theta_applyAT(Ref_E, y);
        else % dense
            Aop.applyA = @(R) theta_applyA_full(E, R);
            Aop.applyAUR = @(R, U) theta_applyA_full(E, R, U);
            Aop.applyAT = @(y) theta_applyAT_full(Ref_E, y);
        end
        
        b = zeros(graph.p, 1);
        
        % construct f
        f = struct();
        f.obj_grad = @theta_obj_grad;
        f.hess = @(R, U, data) 0;
        
        % scaling
        normC = n;
        normA = 1;
        normb = 1;
        normR = 1;
        
        scale.normA = normA; scale.normb = normb;
        scale.normR = normR; scale.normC = normC;
        
        % options
        opts = [];
        opts.verbosity = 1;
        %opts.record_file = ['1.log'];
        opts.record_file = [prob_str, '_', mat2str(test_id), '.log'];
        opts.sub_solver = 1; % default: arnt
        if solver_id == 1
            opts.sub_solver = 3; % 3 for RGBB
        elseif solver_id == 2
            opts.sub_solver = 0; % 0 for adaptive selection
        end
        opts.max_iter = 500;
        opts.tau = 0.8;
        opts.rho = 1.1;
        opts.sigma0  = 100;%sqrt(2*num_E);
        opts.sigma_max = 1e6;
        opts.sigma_min = 1e-2;
        opts.nu_min = 1e-2;
        opts.ALM_step = 1;
        opts.debug = 1;
        opts.comp_eigS = 1;
        if n >= 1500
            opts.use_lobpcg = true;
        else
            opts.use_lobpcg = false;
        end
        opts.adjust_p = 1;
        opts.use_sdpnal = 0;
        opts.gtol_ratio0 = 1e-2;
        opts.gtol0 = 1e-2;
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
        if variant == 0
            opts.subtol_scheme = 1;
            opts.maxitersub = 4;
            if i >= 1 && i <= 2 % for 1dc.128 & 1dc.256, reduce rank(R) non-aggressively
                opts.sigma0 = 10;
                opts.gtol_ratio0 = 1;
                opts.rankR_use_gap = false;
            else
                opts.sigma0 = 100;
                opts.gtol_ratio0 = 1;
            end
            opts.crit = 0.5;
            p = ceil(sqrt(2*graph.p));
            if i == 4 || i == 8 || i == 49 % 1dc.1024 & 1et.1024 & 1dc.2048
                p = ceil(p / 2);
            end
        elseif variant == 1
            p = 2;
            opts.tol = 0.5/n;
            opts.gtol0 = 1;
            opts.sigma0 = 0.1;
            opts.gtol_ratio0 = 1e-2;
            opts.escape_saddle = false;
            opts.comp_eigS = 0;
            opts.crit = 4;
        end
        
        % call ALM_SSN
        [R_,y,S_,Z_,out]=ALM_SSN(n,Aop,b,[],[],'sphere',f,[],p,opts);
        
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
        etaC2 = 0;
        etaC3 = out.deltaC2;
        
        if variant == 0 % standard theta
            
            
            
            
            % table: standalone
        
                time_s = sprintf('%.1f', time);
         
            table_str = [table_str prob_str];
            table_str = [table_str, sprintf(' & %d & %d & %14.5e & %8.1e & %8.1e & %8.1e & %8.1e & %4d / %6.1f / %6.1f & %s & %d & %d', ...
                m, n, pobj, pinf, relgap, etaK2, etaC1, iter, sub_iter,avginit,time_s, X_rank, S_rank)];
            table_str = [table_str ' \\ \hline' newline];
            
            % table: comparison
            
            
         
            
        elseif variant == 1 % rank-two
            % try to recover mss
            alpha1 = round(-pobj);
            [flag, alpha1] = mss_verify(Ref_E, R_, alpha1);
            if flag
                alpha1_s = num2str(alpha1);
            else
                alpha1_s = ['*', num2str(alpha1)];
            end
            
            table_str = [table_str prob_str];
            table_str = [table_str, sprintf(' & %d & %d & %s & %8.1e & %8.1e & %d & %6.1f & %d & %6.1f', ...
                m, n, alpha1_s, pinf, etaC1, iter, time, mss_soln(i,1), mss_soln(i, 2))];
            table_str = [table_str ' \\ \hline' newline];
        end
        
        
    
   
    
    
end

if compress
    rr = 'e([+-])0{0,1}';
    table_str = regexprep(table_str, rr, '$1');
    comp_table_str = regexprep(comp_table_str, rr, '$1');
    detail_table_str = regexprep(detail_table_str, rr, '$1');
end

disp(newline);
disp(table_str);

save_path = strcat(save_root,'test',mat2str(test_type),'_',mat2str(test_id),'.txt');
fid = fopen(save_path,'w+');
fprintf(fid,'%s',table_str);
fclose(fid);

if ~isempty(comp_table_str)
    save_path = strcat(save_root,'comp_test',mat2str(test_type),'_',mat2str(test_id),'.txt');
    fid = fopen(save_path,'w+');
    fprintf(fid,'%s',comp_table_str);
    fclose(fid);
end

if ~isempty(detail_table_str)
    save_path = strcat(save_root,'detail_test',mat2str(test_type),'_',mat2str(test_id),'.txt');
    fid = fopen(save_path,'w+');
    fprintf(fid,'%s',detail_table_str);
    fclose(fid);
end

end
