ns = [512,1024,2048,4096];

num_tiral = 3

for n_idx=1:length(ns)
    for t_idx=1:num_tiral
        n = ns(n_idx);
        u = 1./[1:n];
        u = u'/norm(u);
        V = rand(n,n);
        Sigma = u*u'+2*V'*V;
        save_name = sprintf('../sdp_data/sparse_PCA/random%d_%d.mat',n,t_idx);
        fprintf('%s\n',save_name)
        save(save_name,'Sigma');
    end
end