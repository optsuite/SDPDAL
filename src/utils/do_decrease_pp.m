function [p,R0,M] = do_decrease_pp(R0,p,n,num_decrease,scale_C,Man_handle)


nn = num_decrease; M = Man_handle(p-nn, n, scale_C);
[~,D,V] = svds(R0,p-nn);  
R0 = D*V';


p = p -nn;

end