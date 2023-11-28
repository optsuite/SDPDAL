function ATy_h = nomad_applyAT(y,normA)
n = length(y);
if ~exist('normA', 'var') || isempty(normA)
        normA = 1;
end

ATy_h = @(V) (repmat(V*y,1,n) + sum(V,2)*y')/2/normA ;%- V*diag(diag(y*ones(1,n)));

%V * symm(y *ones(1, length(y))) ;





end





