function [pos, sgn] = cutting_plane_mc(X, K)
%% given X, generate at most K cutting planes that violate
%  X_ij + X_ik + X_jk >= -1
%  X_ij - X_ik - X_jk >= -1
% -X_ij + X_ik - X_jk >= -1
% -X_ij - X_ik + X_jk >= -1

n = size(X, 1);
pos = zeros(4, K);
sgn = zeros(3, K);

ncuts = 0;

for i=1:n
    for j=(i+1):n
        for k=(j+1):n
            v = X(i,j) + X(i,k) + X(j,k);
            if v < -1
                insert(i,j,k,v,1,1,1);
            end
            v = X(i,j) - X(i,k) - X(j,k);
            if v < -1
                insert(i,j,k,v,1,-1,-1);
            end
            v = -X(i,j) + X(i,k) - X(j,k);
            if v < -1
                insert(i,j,k,v,-1,1,-1);
            end
            v = -X(i,j) - X(i,k) + X(j,k);
            if v < -1
                insert(i,j,k,v,-1,-1,1);
            end
        end
    end
end

pos = pos(:,1:ncuts);
sgn = sgn(:,1:ncuts);

    function insert(i, j, k, v, s1, s2, s3)
        if ncuts < K
            ncuts = ncuts + 1;
            pos(:, ncuts) = [i; j; k; v];
            sgn(:, ncuts) = [s1; s2; s3];
        else
            [vmax, imax] = max(pos(4,:));
            if v < vmax
                pos(:,imax) = [i; j; k; v];
                sgn(:,imax) = [s1; s2; s3];
            end
        end
    end
end