function [f, prox_h, prox_h_norm] = pos_L1_obj_prox(R, Z, nuk, data, lbd)

% compute RRZ
if isfield(data, 'RRZ')
    RRZ = data.RRZ;
else
    RRZ = R' * R - Z / nuk;
end

% compute projection of adjoint
W = soft_thres(RRZ, lbd/nuk);
HR = RRZ-W;

f = lbd* sum(abs(W(:)));

% compute prox_h
prox_h = @(V) V * HR;

% compute prox_h_norm
if nargout > 2
    prox_h_norm = norm(HR, 'fro');
end

end