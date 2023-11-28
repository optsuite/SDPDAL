function [F, g, work_out] = obj_maxeig(y, G, delta, work_in)
    n = numel(y);
    % B = G + Ay note: B is sparse
    B = @(x) G * x - y .* x;
    % perform eig
    %work_in.force_exact = 1;
    work_in.degree = 5;
    [V, d, work_out] = eig_inexact(B, work_in, n);
    [ff, gg] = logsumexp_proto(d, delta);
    id = gg > 1e-14;
    
    B = V(:, id) * diag(gg(id)) * V(:, id)';
    
    F = n * ff + sum(y);
    g = -n * diag(B) + 1;
end
