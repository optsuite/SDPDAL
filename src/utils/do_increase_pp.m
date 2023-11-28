function [p,R0,M] = do_increase_pp(S,R0,p,n,num_increase,scale_C,Man_handle,fun_ARNT,varargin)


%tic
%[V,D,flag] = eigs(S,num_increase,'smallestreal','MaxIterations',1e3,'Tolerance',1e-6,'FailureTreatment','keep');
[V,D] = eig(S);
%toc
nn = min(sum(diag(D)<-1e-6), num_increase); 
v = V(:,1:nn); M = Man_handle(p+nn, n, scale_C);
Rt = [R0;zeros(nn,n)];
U = [zeros(p,n);v']; V_step = 1;
R = M.retr(Rt,-U,V_step);

ob0 = feval(fun_ARNT,Rt,varargin{:});
ob1 = feval(fun_ARNT,R,varargin{:});
tt = 1;
while(tt<10 && ob1>1*ob0)
    V_step = V_step * 0.5;
    R = M.retr(Rt,-U,V_step);
    ob1 = feval(fun_ARNT,R,varargin{:});
    tt = tt +1;
end
R0 = R;

p = p + nn;

end