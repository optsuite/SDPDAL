function [blk, At, C, b] = thetaread(fname)
   graph = load(fname);
   E = sortrows(graph.E, 1);
   n = graph.n;
   p = graph.p;

   blk{1,1} = 's';
   blk{1,2} = n;

   C{1} = -sparse(ones(n, n));

   Acell = cell(1, p + 1);

   % tr(X) = 1
   Acell{1} = -speye(n);

   % <E_{ij}, X> = 0, i, j \in E
   for idx = 1:p
       i = E(idx, 1);
       j = E(idx, 2);
       Acell{idx+1} = spconvert([i, j, 1; j, i, 1; n, n, 0]);
   end

   At = svec_sdpnal(blk, Acell, 1);

   b = [-1; zeros(p, 1)];
end
