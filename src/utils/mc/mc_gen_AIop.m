function ret = mc_gen_AIop(Bt)
%% ret = mc_gen_AIop(Bt)
%  Bt - 1-by-nc cell array
%  ret - a struct with
%    .applyA - function handle returns A(R'R) or A(U'R + R'U)
%    .applyAT - function handle returns the operator form of A*(y)

nc = size(Bt, 2);
Bt1 = cell(1, nc);
n = size(Bt{1}, 1);
for i=1:nc
    [ii, jj, vv] = find(Bt{i});
    Bt1{i} = [ii, jj, -vv]';
end
h1 = @(varargin) mex_maxcut_applyA(Bt1, varargin{:});

% concatenate B
Bmat = zeros(3, 3 * nc);
for i=1:nc
    Bmat(:, (3*i-2):3*i) = Bt1{i};
end

h2 = @(y) internal_mc_applyAT(Bmat', n, y);

ret = struct();
ret.applyA = h1;
ret.applyAT = h2;

end

function AITy_h = internal_mc_applyAT(B, n, y)
nc3 = size(B, 1);
nc = nc3 / 3;

Btmp = reshape(B(:, 3), 3, nc);
Btmp = 0.5 * (Btmp .* y');

AITy = spconvert([B(:, [1,2]), reshape(Btmp, nc3, 1); ...
                  B(:, [2,1]), reshape(Btmp, nc3, 1); ...
                  n, n, 0]);

% a sparse-dense matrix production
AITy_h = @(V) V * AITy;

end