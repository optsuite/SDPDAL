function flag = etaK2_stagnate_check(list, itr, stp)

if itr <= stp + 2
    flag = false;
    return;
end

ratio = list(itr-stp-1:itr-1) ./ list(itr-stp-2:itr-2);

if sum(ratio > 1) > 0.5 * stp && min(ratio) > 0.97 && max(ratio < 1.5)
    flag = true;
else
    flag = false;
end