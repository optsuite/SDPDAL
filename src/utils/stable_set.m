function [stablesetsize, mss] = stable_set(E, R)
AX = theta_applyA(E, R);
u = R(1, :);
v = R(2, :);

eps = max(abs(AX));
eps1 = 0.5*eps; eps2 = sqrt(eps1);
eps3 = sqrt(eps*(1.0 + 0.25*eps)); eps4 = sqrt(2*eps);

set1 = u.*v < 0 & abs(u) > eps2 & abs(v) > eps2; size1 = nnz(set1);
set2 = u.*v > 0 & abs(u) > eps2 & abs(v) > eps2; size2 = nnz(set2);
set3 = abs(v) < eps1 & abs(u) >  eps3 & abs(u) + abs(u) > eps4; size3 = nnz(set3);
set4 = abs(u) < eps1 & abs(v) >  eps3 & abs(v) + abs(v) > eps4; size4 = nnz(set4);
%if(size1*size2 < 0)
%    tempervec = set1 + set2 + set3 + set4;
%end

if(size1 > size2)
    if(size3 > size4)
        stablesetsize = size1 + size3;  setw = set1 + set3;
    else
        stablesetsize = size1 + size4;  setw = set1 + set4;
    end
else
    if(size3 > size4)
        stablesetsize = size2 + size3;  setw = set2 + set3;
    else
        stablesetsize = size2 + size4;  setw = set2 + set4;
    end
end

% verify setw
if any( setw(E(:,1)) .* setw(E(:,2)) == 1 )
    mss = [];
    stablesetsize = 0;
else
    mss = setw;
end

end
