function [u, BTu] = compute_u(manifold, gradR, R, normR)

if strcmp(manifold, 'sphere')
    n = size(R,2);
    u = sum(sum(gradR.*R))/2/normR^2;
    BTu = eye(n)*u;
end

if strcmp(manifold, 'oblique') 
    u = diag(R'*gradR)/2/normR^2;
    BTu = diag(u);
end

if strcmp(manifold, 'oblique_shape') 
    n = sqrt(size(R,2));
    RgradR = zeros(n,1);
    for i=1:n
          RgradR(i) = sum(sum( R(:,(i-1)*n+1:i*n).* gradR(:,(i-1)*n+1:i*n))) ;
    end
    u = RgradR/2/normR^2;
    BTu = zeros(n^2,n^2);
    for i=1:n
        BTu((i-1)*n+1:i*n, (i-1)*n+1:i*n) = u(i)*eye(n);
    end
end

if strcmp(manifold, 'stiefel') 
    n = sqrt(size(R,2));
    RgradR = zeros(n,n);
    for i=1:n
        for j=1:n
          RgradR(i,j) = sum(sum( R(:,(i-1)*n+1:i*n).* gradR(:,(j-1)*n+1:j*n))) ;
        end
    end
    u = (RgradR + RgradR')/4/normR^2;
    BTu = zeros(n^2,n^2);
    for i=1:n
        for j=1:n
        BTu((i-1)*n+1:i*n, (j-1)*n+1:j*n) = u(i,j)*eye(n);
        end
    end
    
end

