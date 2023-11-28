function Test_rcp(test_type, test_id, rank_thres, compress)
%clear;
% min <C,X> s.t. Xe=e, tr(X) = K, X>=0, X \succeq 0 


if nargin<1; test_type =1; end
if nargin<2; test_id = 1;   end
if nargin<3; rank_thres = 1e-5;   end
if nargin<4; compress = true; end
% 0: full 1-4: small problems
% 5-14: big problems
% -1: debug
addpath ../
file_dir = '..';
root_dir = '../..';
addpath(file_dir);
addpath(genpath([root_dir '/SDPNAL+v1.0']));



addpath(genpath(strcat(file_dir,'/src')));
addpath(strcat(file_dir,'/ARNT/solvers'));
addpath(strcat(file_dir,'/ARNT/manifolds-manopt'));
data_dir = [file_dir,'/sdp_data/'];



if test_type ==0
    [files, nlist] = rcpprobs;
elseif test_type>=1 && test_type<=4 && rem(test_type,1)==0
    [files, nlist] = rcpprobs_small(test_type);
elseif test_type>=5 && test_type<=14 && rem(test_type,1)==0
    [files, nlist] = rcpprobs_big;
elseif test_type == 40
    [files, nlist] = rcpprobs_200;
elseif test_type == 39
    [files, nlist] = rcpprobs_segment;
elseif test_type == -1
    [files, nlist] = rcpprobs_debug;
else
    files = [];
    fprintf('wrong test_type.\n');
end
file_len = length(files);

prob_range = 1:file_len;

table_str = '';
detail_table_str = '';

% profiling
etime_ALM = [];
etime_NAL = [];

save_root = strcat(file_dir,'/results/rcp/');
if ~exist(save_root,'dir')
    mkdir(save_root)
end

save_root_mat = strcat(save_root,'mat/');
if ~exist(save_root_mat,'dir')
    mkdir(save_root_mat)
end

save_root_eig = strcat(save_root,'eig/');
if ~exist(save_root_eig,'dir')
    mkdir(save_root_eig)
end

save_root_res = strcat(save_root,'res/');
if ~exist(save_root_res,'dir')
    mkdir(save_root_res)
end
%prob_range = [2];
for i = prob_range
    relpath = files{i};
    [~, basename1, basename2] = fileparts(relpath);
    %basename = [basename1, basename2];
    file = [data_dir, 'rcp/', relpath, '.data'];
    for size_idx = 1:length(nlist{i})
        [n0, postfix] = deal(nlist{i}{size_idx}{:});
        if test_type<=40
            etime_ALM_sub = zeros(1, 10);
            etime_NAL_sub = zeros(1, 10);
            for K = 2:11
                basename = [basename1, basename2, postfix, '.', num2str(K)];
                
                W0 = importdata(file);
                if isstruct(W0)
                    W0 = W0.data;
                end
                
                X = W0(1:n0,:);
                m = n0;
                
                scale_C = 1e0;
                p = 1*min(ceil(sqrt(2 * m)), 200);
                p = max(2*K,min(p,4*K))*2;
                
                C = -1*(X*X')/1;
                
                % scale parameter
                normC = (norm(C,'fro'));
                normb = sqrt(K+n0);
                normA = sqrt(n0*(n0+2)/2);
                normR = sqrt(K*normA/normb);
                
                
                % scale parameter version2
                % Ai has norm <= 1
                % tr(X) = 1
                % C has norm <= 1
                normC = max(1, norm(C, 'fro'));
                normA = sqrt((n0+1)/2);%sqrt((n0+2))/2;
                normb = K*normA;
                normR = 1;
                
                % Original version
                %                 normC = (norm(C,'fro'));
                %                 normb = K;
                %                 normA = 1;
                %                 normR = sqrt(K*normA/normb);
                
                Y = X/sqrt(normC/1);
                
                
                
                scale.normA = normA; scale.normb = normb;
                scale.normR = normR; scale.normC = normC;
                
                
                
                % construct Aop
                Aop.applyA = @(R) nomad_applyA(R,normA);   % A(R^TR)
                Aop.applyAUR = @(R,U) nomad_applyA(R,normA,U); % A(R^TU + U^TR)
                Aop.applyAT = @(y) nomad_applyAT(y,normA);  % A^T(y)
                
                
                
                rng(2333);
                b = ones(m,1)/normb;
                
                % construct f
                
                f = struct();
                f.obj_grad = @nomad_obj_grad;
                f.hess = @(R, U, data, C) 0;
                f.data = {Y};
                
                
                lbd = 0;
                
                % construct h
                h = struct();
                h.obj_prox = @pos_box_obj_prox;
                h.hess = @pos_box_hess;
                h.data = {lbd};
                
                gtol_ratio0 = 1e-3 / K;
                opts = [];
                opts.verbosity = 2;
                %opts.record_file = [basename '.log'];
                if test_id == 1
                    opts.sub_solver = 3; % 3 for RGBB
                elseif test_id == 2
                    opts.sub_solver = 0; % 0 for adaptive selection
                end
                opts.max_iter = 500;
                opts.tau = 0.8;
                opts.rho = 1.1;
                %opts.sigma0 = norm(C,'fro')/100;
                % tune nu0_fac
                if strcmp(basename1, 'housing')
                    nu0_fac = 5;
                elseif strcmp(basename1, 'abalone')
                    nu0_fac = 5;
                elseif strcmp(basename1, 'segment')
                    if strcmp(postfix, '-small')
                        nu0_fac = 100;
                        if K == 5 || K == 6
                            nu0_fac = 400; 
                        end
                    elseif strcmp(postfix, '-medium')
                        nu0_fac = 200;
                    elseif strcmp(postfix, '-large')
                        nu0_fac = 200;
                    end
                    % smaller gtol_ratio0 for escaping saddle
                    if K >= 6; gtol_ratio0 = 1e-4 / K; end
                elseif strcmp(basename1, 'spambase')
                    if strcmp(postfix, '-small')
                        nu0_fac = 50;
                    elseif strcmp(postfix, '-medium')
                        nu0_fac = 150;
                    elseif strcmp(postfix, '-large')
                        nu0_fac = 250;
                      
                    end
                elseif strcmp(basename1, 'soybean-large')
                    nu0_fac = 50;
                else
                    nu0_fac = 5;
                    %opts.sub_solver = 3;
                end
                opts.nu0 = nu0_fac/K;
                opts.sigma0  = nu0_fac/K;
                opts.sigma_min = 1e-2;
                opts.nu_min = 1e-2;
                opts.nu_max = opts.nu0 * 1000;
                opts.use_mex = 1;
                opts.ALM_step = 1;
                opts.debug = 1;
                opts.comp_eigS = 1;
                opts.adjust_p = 1;
                opts.use_sdpnal = 0;
                opts.gtol_ratio0 = gtol_ratio0;
                %opts.gtol_ratio0 = 2;
                opts.flag_decrease_p = true;
                opts.favor_deltak = true;
                opts.disable_auto_sigma0 = 1;
                opts.scale_C = scale_C;
                opts.scale = scale;
                opts.tol = 1e-6;
                opts.escape_saddle = true;
                %opts.early_stop = true;
                opts.crit = 1;
                if opts.sub_solver ~= 1
                    opts.verbosity = 1;
                end
                
                % construct initial point
                ss = 5;   % better than 1
                M = (ss*p+1-ss*K)/(K*(m-1)*(ss*p+1));
                X1 = M*ones(p,m);
                for k=1:p
                    X1(k,k) =  1*ss/(ss*p+1);
                end
                X2 = X1(1:p,p+1:m)';
                N = 1/(ss*p+1)/(m-p);
                O = (1/K - M*p - N)/(m-p-1);
                X3 = O*ones(m-p,m-p);
                for k=1:m-p
                    X3(k,k) =  N;
                end
                Q = [X1;X2,X3];
                [R,D] = eigs(Q,p);
                R = R*sqrt(D);
                % opts.R0 = R'/norm(R,'fro');
                opts.R0 = R';
                
                % initial point (version 2)
                idx = randi(p, 1, n0);
                R1 = zeros(p, n0);
                for ii = 1:n0
                    R1(idx(ii), ii) = 1;
                end
                Rs = sum(R1, 2);
                R1 = R1 ./ sqrt(Rs * p);
                opts.R0 = R1;
                
                [R_,y,S_,Z_,out]=SDPDAL(m,Aop,b,[],[],'sphere',f,h,p,opts);
                
                
                
                R_ = sqrt(1)*R_;
                X = {};
                X{1} = (R_'*R_);
                S = {};
                S{1} = S_;
                Z = {};
                Z{1} = Z_;
                %                 [obj,X,~,y,S,Z,~,~,out,runhist] = ...
                %                     sdpnalplus(blk,At,C,b,L,U,[],[],[],OPTIONS);
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
               
                
                
                % X_rank = sum(X_eig-max(X_eig)*rank_thres>0);
                % S_rank = sum(S_eig-max(S_eig)*rank_thres>0);
                %                 X_rank = sum(X_eig-max(X_eig)*relgap>0);
                %                 S_rank = sum(S_eig-max(S_eig)*relgap>0);
                
                X_rank = sub_comp_rank(X_eig,struct('thres',rank_thres));
                S_rank = sub_comp_rank(S_eig,struct('thres',rank_thres));
                
                etaK2 = out.etaK2;
                etaC1 = out.deltaC; avginit = out.HessianVector;
                sub_iter = out.sub_iter;
                etaK1 = 0;
                etaC2 = abs(sum(sum(X_item .* Z_))) / (1 + norm(X_item, 'fro') + norm(Z_, 'fro'));
                %sdprcp = load([save_root, 'res_sdpnal/', basename,'_res_',mat2str(1),'.mat']);
                p = out.p;
                marker = '';
                if out.flag >= 2; marker = '*'; end
                
                
                prob_str1 = basename1(1:2);
                if strcmp(postfix, '')
                    prob_str2 = '';
                else
                    prob_str2 = postfix(1:2);
                    prob_str2(1) = '.';
                end
                prob_str = [prob_str1, prob_str2, '.', num2str(K)];
                table_str = [table_str marker prob_str];
                table_str = [table_str, sprintf('& %d & %.2e & %.1e & %.1e & %.1e & %.1e & %.1e & %d / %.1f / %.1f & %.1f & %d & %d', ...
                    n, pobj, pinf, relgap, etaK2, etaC1, etaC2, iter, sub_iter,avginit,time, X_rank, S_rank)];
                table_str = [table_str ' \\ \hline' newline];
                
                % print the detail of sdpnal
                
                
                
                etime_ALM_sub(K-1) = time;
                %etime_NAL_sub(K-1) = sdprcp.time;
                
            end
            etime_ALM = [etime_ALM, etime_ALM_sub];
            etime_NAL = [etime_NAL, etime_NAL_sub];
            
            
        end
    end
    
    if compress
        rr = 'e([+-])0{0,1}';
        table_str = regexprep(table_str, rr, '$1');
        
    end
    
    disp(newline);
    disp(table_str);
    
   
    
    
    
    
    
end
