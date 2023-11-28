function [f, prox_h, prox_h_norm] = pos_non_obj_prox(R, Z, nuk, data,lb)



% after projection, f is always 0
f = 0;

% compute prox_h
prox_h = @(V) V * 0;

% compute prox_h_norm
if nargout > 2
    prox_h_norm = 0;
end

end