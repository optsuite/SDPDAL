function [TX] = compute_TX_box(R, R0, nuk,i,lb)
% compute  the i-th column of TX.
%  TX = R'*R + Z/nuk - prox_h(R'*R + Z/nuk), where Z = R0'*R0 -prox_h(R0'*R0);



rrz = R'*R(:,i) +  min(R0'*R0(:,i),lb)/nuk;
TX = min(rrz,lb);