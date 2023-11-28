function y = nomad_applyA(R,normA,U)
 
if ~exist('normA', 'var') || isempty(normA)
        normA = 1;
end
    
if nargin<3
    y = R' * (sum(R,2))/normA;
else
    y =(R' * (sum(U,2)) + U' * (sum(R,2)))/normA;
end
