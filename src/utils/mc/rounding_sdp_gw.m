function [cut, x] = rounding_sdp_gw(C, R, ntimes)

cut_best = 0;
p = size(R, 1);
for i = 1:ntimes
    r = randn(p, 1);
    x = sign(R' * r);
    cut = x' * (C * x);
    if i == 1 || cut > cut_best
        cut_best = cut;
        x_best = x;
    end
end

cut = cut_best;
x = x_best;

end