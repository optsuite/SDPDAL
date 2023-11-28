function [cut, x] = rounding_local_1_opt(C, x)

n = size(x, 1);

i = 1;
while i <= n 
    cut = x' * (C * x);
    x1 = x;
    % perform flipping
    x1(i) = -x1(i);
    cut1 = x1' * (C * x1);
    if cut1 > cut
        x = x1;
        cut = cut1;
        % restart from i = 1
        i = 1;
    else
        i = i + 1;
    end
end

end