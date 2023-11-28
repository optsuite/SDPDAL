function Test_ncm()
%% Find nearest correlation matrix
%  \min 0.5*\|H.*(X-G)\|_F^2, s.t. diag(X) = e, X>L, X is psd.


rank_thres = 1e-5;   
% 0: full 1: small problems
% 2-30: big problems
% -1: debug



file_dir = '..';
addpath(file_dir);



addpath(genpath(strcat(file_dir,'/src')));
addpath(strcat(file_dir,'/ARNT/solvers'));
addpath(strcat(file_dir,'/ARNT/manifolds-manopt'));
data_dir = [file_dir '/sdp_data/ncm/'];





table_str = '';



fileFolder=fullfile(data_dir);
dirOutput=dir(fullfile(fileFolder,'*.mat'));
fileNames={dirOutput.name};


k = length(fileNames);


load([file_dir,'/sdp_data/divergentGH.mat'],'H');
tmp0 = ones(110,110); d  = 10;
H0 = kron(tmp0,H);





basename = 'Leukemia';
alpha = 0.1; lbd = -0.3;
load([data_dir,basename],'S');
n = length(S);
%                 H = H0([1:n],[1:n]);   choose H
%                 H = (H+H')/2;
H = 1;  % without H
G = S([1:n],[1:n]);


clear S;
rng(2223);


R = rand(d,n); F = rand(d,n);E = R'*F;
G = (1-alpha)*G + 1*alpha*E;
G = 0.5*(G+G');

G = G - diag(diag(G)) + eye(n);
tempH = H;
Hscale = sqrt(norm(tempH.*tempH.*G,'fro'));
H = H/(Hscale);


scale_C = 1;
% scale parameter version2
% Ai has norm <= 1
% diag(X) = 1
% C has norm <= 1
normC = 1;
normA = 1;
normb = 1;
normR = sqrt(normA/normb);





scale.normA = normA; scale.normb = normb;
scale.normR = normR; scale.normC = normC;



% construct Aop
Aop.applyA = @(R) 0; 
Aop.applyAUR = @(R,U) 0;
Aop.applyAT = @AT;
b = 0;



rng(2223);


% construct f

f = struct();
f.obj_grad = @ncm_obj_grad;
f.hess = @ncm_hess;
f.data = {H,G};
f.conjugate = @ncm_conjugate;




% construct h
h = struct();

h.obj_prox = @pos_box_obj_prox;
h.hess = @pos_box_hess;

h.data = {lbd};
h.conjugate = @(Z,lbd) lbd*sum(sum(Z));

opts = [];
opts.verbosity = 1;


opts.sub_solver = 1;  % 0:ÇÐ»»£» 1£ºarnt 3£ºgradient
opts.max_iter = 100;
opts.tau = 0.8;
opts.rho = 1.2;
%opts.sigma0 = norm(C,'fro')/100;
% tune nu0_fac

opts.nu0 = 0.1;
opts.sigma0  = 1;
opts.sigma_min = 1e-4;
opts.nu_min = 1e-4;
%opts.nu_max = 20;
opts.use_mex = 1;
opts.ALM_step = 1;
opts.debug = 1;
opts.comp_eigS = 1;
opts.adjust_p = 1;
opts.use_sdpnal = 0;
opts.gtol_ratio0 = scale_C/1e3;
%opts.gtol_ratio0 = 2;
opts.flag_decrease_p = false;
opts.favor_deltak = true;
opts.disable_auto_sigma0 = 1;
opts.scale_C = scale_C;
opts.scale = scale;
opts.tol = 1e-6;
opts.escape_saddle = 0  ;
opts.increase_p = 1;
opts.decrease_p = 1;
opts.early_stop = true;
opts.gtol_min = 1e-6;
opts.gtol_decrease = 0.6;
opts.gtol0 = 1e-1;


p = 100;
[V,D] = eigs(G,p,'largestreal');
D(D<0) = 0;
R = (D.^0.5)*V';

for kk = 1:n
    R(:,kk) = R(:,kk)/norm(R(:,kk))*normR;
end
opts.R0 = R;

[R_,y,S_,Z_,out]=SDPDAL(n,Aop,b,[],[],'oblique',f,h,p,opts);



X = {};
X{1} = (R_'*R_);
S = {};
S{1} = S_;
Z = {};
Z{1} = Z_;
XZ = X{1}-Z{1}; XZ(XZ<lbd)=lbd;


iter = out.iter;
time = out.time;
pinf = out.deltak;
dinf = 0;
pobj = 0.5*norm(H.*(X{1} - G), 'fro')^2 * Hscale^2;% - mu/2*norm(H.*X{1},'fro')^2;%out.pobj;
dobj = out.dobj;
relgap = out.relgap;

X_item = X{1}; %X_item(X_item<lbd)=lbd;
[n,~] = size(X_item);




X_eig = eig(X_item);
S_eig = eig(S{1});


S_eig_neg = S_eig.*(S_eig<0);
X_eig_neg = X_eig.*(X_eig<0);

X_rank = sub_comp_rank(X_eig,struct('thres',rank_thres));
S_rank = sub_comp_rank(S_eig,struct('thres',rank_thres));

etaK2 = out.etaK2;
etaC1 = out.deltaC; avginit = out.HessianVector;
sub_iter = out.sub_iter;
etaK1 = 0;
etaC2 = norm(X{1}-XZ,'fro')/(1+norm(X{1},'fro')+norm(Z{1},'fro'));
nrmG = out.nrmG;

p = out.p;

table_str = [table_str basename '-' mat2str(alpha),'-',mat2str(lbd)];
table_str = [table_str, sprintf('& %d & %d & %.2e & %.2e & %.2e & %.2e & %.2e & %d & %.1f & %.1f & %.1f & %d & %d', ...
    p,n, pobj, relgap, etaK2, etaC1, etaC2, iter, sub_iter,avginit,time, X_rank, S_rank)];
table_str = [table_str ' \\ \hline' newline];








disp(newline);
disp(table_str);



    function re = AT(y)
        re = @(V) 0*y;
        
    end

    function re = ncm_hess(R, U, data, H,G)
        if isfield(data,'UtR')
            UtR = data.UtR;
        else
            UtR = R'*U; UtR = UtR' + UtR;
        end
        if(~isempty(H))
            re = R * ((H.^2).*UtR);
        else
            re = R * UtR;
        end
    end




end
