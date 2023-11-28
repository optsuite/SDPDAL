function [TF, true_alpha] = rss_verify(Ref_E, R, a)
%% verify a rank-two solution to be valid

nrmR = sqrt(sum(R .* R, 1));
[~, id] = sort(nrmR, 'descend');
for i = a:-1:1
    E_sub = Ref_E(id(1:i), id(1:i));
    if nnz(E_sub) == 0
        break;
    end
end

% check i
true_alpha = i;
TF = (i == a);
end
