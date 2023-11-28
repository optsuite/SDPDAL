function mat2sparse(A, filename)

fid = fopen(filename, 'w');
if issymmetric(A)
    % only save upper for symmetric A
    A = triu(A);
end

[ii, jj, vv] = find(A);

fprintf(fid, '%d %d\n', size(A,1), nnz(A));
% write data
fprintf(fid, '%d %d %g\n', [ii, jj, vv]');

fclose(fid);