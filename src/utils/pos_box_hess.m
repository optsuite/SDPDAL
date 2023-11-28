function h = pos_box_hess(R, Z, nuk, U, data,lb)

% compute RRZ
if isfield(data, 'RRZ')
    RRZ = data.RRZ;
else
    RRZ = R' * R - Z / nuk;
end
if isfield(data,'UtR')
    UtR = data.UtR;
else
    UtR = mexfun_mex_almssn_compute_pgU(U,R);
end
%UtR = U'*R; UtR = UtR + UtR';
pgU = UtR.* (RRZ < lb);

%pgU = pgU+pgU';
%pgU = (UtR + UtR') .* a;

% pgU = zeros(200);
% index = (RRZ < lb);
% for i=1:200
%     for j=1:200
%         if(index(i,j)==1)
%             pgU(i,j) = U(:,i)'*R(:,j) + R(:,i)'*U(:,j);
%         end
%     end
% end


% compute H2
%h = mexfun_mex_almssn_compute_RpgU(R,pgU) ;
h = R*pgU;
end