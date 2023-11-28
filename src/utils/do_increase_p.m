function [R_new, M,p, outinfo] = do_increase_p(R, S, fun, varargin)
	[p, n] = size(R);
	[Vmin,~,flag] = eigs(S,1,'smallestreal','MaxIterations',1e3,'Tolerance',1e-6,'FailureTreatment','keep');
	if flag
		fprintf('\n Eigenvector not converge.');
	end
	p = p+1;
	R = [R;zeros(1,n)];
	M = spherefactory(p,n);
	V = [zeros(p-1,n);Vmin'];
	V_step = 1/(p-1);
	eta_ = 0.2;
	mu = 1e-4;
	uSu = Vmin'*S*Vmin;
	f_old = feval(fun,R,varargin{:});
	nls = 1;
	while 1
		R_new = M.retr(R,V,V_step);
		f_new = feval(fun,R_new,varargin{:});
		if f_new-f_old<mu*uSu*V_step^2 || nls>=10
			break;
		end
		V_step = V_step*eta_;
		nls = nls+1;
	end
	outinfo.nls = nls;

end