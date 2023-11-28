function Hg = getBg_Lbfgs(g,S,Y)
% This function returns the approximate inverse Hessian multiplied by the gradient, H*g
% Input
%   S:    Memory matrix (n by k) , s{i}=x{i+1}-x{i}
%   Y:    Memory matrix (n by k) , df{i}=df{i+1}-df{i}
%   g:    gradient (n by 1)
%   hdiag value of initial Hessian diagonal elements (scalar)
% Output
%   Hg    the the approximate inverse Hessian multiplied by the gradient g
% Notice
% This funcion getHg_lbfgs is called by LBFGS_opt.m.
% Ref
%   Nocedal, J. (1980). "Updating Quasi-Newton Matrices with Limited Storage".
%   Wiki http://en.wikipedia.org/wiki/Limited-memory_BFGS
%   two loop recursion

    k = length(S); 
    [r,n] = size(g);
    SS = zeros(k,k);
    L = zeros(k,k);
    D = zeros(k,k);
    hdiag = sum(sum(S{end}.*Y{end}))/sum(sum(S{end}.*S{end}));
    for i = 1:k
        for j = i:k
            SS(i,j) = hdiag*sum(sum(S{i}.*S{j}));
            SS(j,i) = SS(i,j);
            if(i>j)
                L(i,j) = sum(sum(S{i}.*Y{j}));
            else
                D(i,i) = sum(sum(S{i}.*Y{j}));
            end
        end
    end
    
    STg = zeros(k,1); YTg = zeros(k,1);
    for i = 1:k
        STg(i) = hdiag*sum(sum(S{i}.*g));
        YTg(i) = sum(sum(Y{i}.*g));
    end
    
    C = [SS,L;L',-D];
    
    Cg = C\[STg;YTg];
    
    SYg = zeros(r,n);
    for i=1:k
        SYg = SYg + hdiag*Cg(i)*S{i} + Cg(k+i)*Y{i};
    end
    
    Hg = hdiag*g - SYg;

end % end of getHg_lbfgs