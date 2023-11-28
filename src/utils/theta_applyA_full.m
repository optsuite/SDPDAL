function AX = theta_applyA_full(E, varargin)

% for R'R
if nargin < 3
    X = varargin{1}' * varargin{1};
else % for R'U + U'R
    X = varargin{1}' * varargin{2};
    X = X + X';
end

AX = mex_theta_applyA_full(E, X);

end
