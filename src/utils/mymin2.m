  function [fmin, I,J] = mymin2(A)
      % A = A + 10*speye(size(A,1));
       [fmin_row, index] = min(A,[],2);
       [fmin, index2] = min(fmin_row);
       I = index2;
       J = index(I);
  end