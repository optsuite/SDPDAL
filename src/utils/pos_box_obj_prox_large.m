function [f, UTX, prox_h_norm] = pos_box_obj_prox_large(R, R0, U,nuk, data,lb)



%% Compute the Moreau envelope of box constraint.
%   f = h(prox_h(R'*R + Z/nuk));
%   TX = R'*R + Z/nuk - prox_h(R'*R + Z/nuk), where Z = R0'*R0 -prox_h(R0'*R0);
%   UTX = U * TX; 
%   prox_h_norm = \|TX\|_F



RRZ_minus_vector = zeros(size(R,2),1);


if isfield(data, 'UTX')
    UTX = data.UTX;
else
    UTX = zeros(size(R));
    for i=1:size(R,2)
        %compute the i-th column of TX
        rrz_minus = compute_TX_box(R, R0, nuk,i,lb);
        if nargout > 2
            RRZ_minus_vector(i) = norm(rrz_minus);
        end
        
        UTX(:,i) = U*rrz_minus;
    end
end



% compute projection of adjoint

% after projection, f is always 0
f = 0;

% compute prox_h


% compute prox_h_norm
if nargout > 2
    prox_h_norm = norm(RRZ_minus_vector);
end

end