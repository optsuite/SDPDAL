function [f, prox_h, prox_h_norm] = pos_box_obj_prox(R, Z, nuk, data,lb)

% compute RRZ
if isfield(data, 'RRZ')
    RRZ = data.RRZ;
elseif isfield(data,'RR')
        RRZ = data.RR - Z / nuk;
    else
        RRZ = R' * R - Z / nuk;
end

% compute projection of adjoint

RRZ_minus = (min(RRZ-lb, 0));

% after projection, f is always 0
f = 0;

% compute prox_h
prox_h = @(V) V * RRZ_minus;

% compute prox_h_norm
if nargout > 2
    prox_h_norm = norm(RRZ_minus, 'fro');
end

end