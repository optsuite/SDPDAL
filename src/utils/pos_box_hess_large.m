function h = pos_box_hess_large(R, R0, nuk, U, data,lb)

% TX = R'*R + Z/nuk - prox_h(R'*R + Z/nuk), where Z = R0'*R0 -prox_h(R0'*R0);
% h = R * \partial TX * [R'*U + U'*R];
% 
h = zeros(size(R));
for i = 1:size(R,2)
    TXi = compute_TX_box(R, R0, nuk, i, lb);
    pgU = (R'*U(:,i) + U'*R(:,i)).*(TXi < lb);
    h(:,i) = R*pgU;
end



end