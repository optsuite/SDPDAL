function Y= soft_thres(X,lbd)
% soft threshold
    Y = sign(X).*max(abs(X)-lbd,0);
end