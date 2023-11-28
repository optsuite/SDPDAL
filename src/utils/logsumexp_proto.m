function [f, g, G] = logsumexp_proto(x, delta)
    a = max(x);
    exp_x = exp((x - a) / delta);
    % function value
    f = a + delta * log(sum(exp_x));
    
    % gradient
    if nargout >= 2
        g = exp_x / sum(exp_x);
    end
    
    % hessian
    if nargout >= 3
        % TODO
        G = 0;
    end
end