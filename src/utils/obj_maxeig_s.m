function [F, g, work_out] = obj_maxeig_s(y, G, work_in)
    n = numel(y);
    alpha = n;
    % B = G + Ay note: B is sparse
    B = @(x) y .* x - G * x;
    % perform eig
    work_in.degree = 5;
    work_in.itr = 1;
    %work_in.force_exact = 1;
    [V, d, work_out] = eig_inexact(B, work_in, n);
    work_out.d = d;
    ff = max(d);
    id = d > ff - 1e-5;
    [~, s] = logsumexp_proto(d(id), 1e-5);
    S = spdiags(s, 0, sum(id), sum(id));
    
    B = V(:, id) * S * V(:, id)';
    
    F = alpha * max(ff, 0) - sum(y);
    
    if ff > 0
        g = alpha * diag(B) - 1;
    else
        g = -ones(n, 1);
    end
end