function obj = ncm_conjugate(X, H, G)


if(~isempty(H))
    obj = 0.5*norm(X./H,'fro')^2 + sum(sum(G.*X));
else
    obj = 0.5*norm(X,'fro')^2 + sum(sum(G.*X));
end


% if(isempty(H))
%     H = ones(size(G));
% end
% Z = (X + (H.^2).*G)./((1-mu)*H.^2);
% obj = sum(sum(X.*Z)) - 0.5*norm(H.*(Z-G),'fro')^2 + 0.5*mu*norm(H.*Z,'fro')^2;



end