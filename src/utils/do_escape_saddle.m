function [R_new, flag] = do_escape_saddle(Seig, Q1, R1,rankR,M, fun, varargin)
    [p, n] = size(R1);
    if ~isstruct(Seig)
        d = ceil((p - rankR)/2);  [Vmin,uSu,~] = eigs(Seig,d,'smallestreal'); 
        Seig = struct(); Seig.d = diag(uSu); Seig.V = Vmin;
    end
    uSu = min(0, Seig.d);

    % how many directions can be used for escaping saddle?
    d = min(ceil((p - rankR)/2), sum(uSu < 0));

    % if there is no descent direction or S is not available
    if d == 0 || isempty(Seig)
        R_new = 0;
        flag = 1;
        return;
    end
        
    uSu = -norm(uSu,'fro');    
    Vmin = Seig.V(:,1:d);
    flag = 1;
	%R1(end,:) = 0;
	V = zeros(p,n);
	V(end-d+1:end,:) = -Vmin';
	V_step = 1/rankR;
	eta_ = 0.5;
	mu = 1e-3;
	%uSu = Vmin'*S*Vmin;
	f_old = feval(fun,R1,varargin{:});
	nls = 1;

	while 1
		R_new = M.retr(R1,V,V_step);
		f_new = feval(fun,R_new,varargin{:});
		if f_new-f_old<mu*uSu*V_step^2 || nls>=10
            flag = 0;
			break;
		end
		V_step = V_step*eta_;
		nls = nls+1;
	end

	R_new = Q1*R_new;

end
