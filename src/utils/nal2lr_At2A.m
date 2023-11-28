function [A] =  nal2lr_At2A(At)
	[n_,m] = size(At{1});
	n = 0.5*(sqrt(8*n_+1)-1);
	[indx,indy,entryz] = find(At{1});
	indx_row = ceil(0.5*(sqrt(8*indx+1)-1));
	indx_col = indx-indx_row.*(indx_row-1)/2;

	assert(sum(indx_row-n>0)==0)
	assert(sum(indx_col-n>0)==0)
	A = sparse(m,n^2);

	for j = 1:length(indx)
		y = indy(j);
		row = indx_row(j);
		col = indx_col(j);
		z = entryz(j);
		if row==col
			A(y,(row-1)*n+col) = z;
		else
			A(y,(row-1)*n+col) = z/sqrt(2);
			A(y,(col-1)*n+row) = z/sqrt(2);
		end
	end
end