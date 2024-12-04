function res = orth_basis_legendre(x, deg, drv)
%ORTH_BASIS_LEGENDRE Compute values or derivatives of Legendre basis.
%   This function calculates the Legendre basis polynomials (orthogonal
%   over [-1/2, 1/2]) and their derivatives up to degree 4.
%
%   Syntax:
%      res = ORTH_BASIS_LEGENDRE(x, deg)
%      res = ORTH_BASIS_LEGENDRE(x, deg, drv)
%
%   Inputs:
%      x   - Scalar or array, points within the interval [-1/2, 1/2].
%      deg - Degree of the Legendre polynomial (0 <= deg <= 4).
%      drv - (Optional) Order of the derivative (0 or 1). Default is 0.
%
%   Output:
%      res - Values of the Legendre polynomial or its derivative at `x`.
%
%   Example:
%      % Compute the degree-3 Legendre polynomial at x = 0.1
%      res = orth_basis_legendre(0.1, 3);
%
%      % Compute the derivative of the degree-3 polynomial at x = 0.1
%      res = orth_basis_legendre(0.1, 3, 1);
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

if nargin < 3
    drv = 0; % Default: no derivative
end

% Validate inputs
if deg < 0 || deg > 4

end
if drv < 0 || drv > 1
    error('Derivative order must be 0 or 1.');
end

% Dispatch to helper function
res = eval_basis(x, deg, drv);
end

function res = eval_basis(x, deg, drv)
% Dispatch to degree-specific function
switch deg
    case 0
        res = eval_deg_0(x, drv);
    case 1
        res = eval_deg_1(x, drv);
    case 2
        res = eval_deg_2(x, drv);
    case 3
        res = eval_deg_3(x, drv);
    case 4
        res = eval_deg_4(x, drv);
end
end

function res = eval_deg_0(x, drv)
switch drv
    case 0
        res = ones(size(x));
    case 1
        res = zeros(size(x));
end
end

function res = eval_deg_1(x, drv)
switch drv
    case 0
        res = x;
    case 1
        res = ones(size(x));
end
end

function res = eval_deg_2(x, drv)
switch drv
    case 0
        res = x.^2 - 1 / 12;
    case 1
        res = 2 * x;
end
end

function res = eval_deg_3(x, drv)
switch drv
    case 0
        res = (x.^2 - 3 / 20) .* x;
    case 1
        res = 3 * x.^2 - 3 / 20;
end
end

function res = eval_deg_4(x, drv)
switch drv
    case 0
        res = (x.^2 - 3 / 14) .* x.^2 + 3 / 560;
    case 1
        res = 4 * x .* (x.^2 - 3 / 14);
end
end