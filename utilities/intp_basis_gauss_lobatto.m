function [res] = intp_basis_gauss_lobatto(x, idx, deg, drv)
%INTP_BASIS_GAUSS_LOBATTO Compute interpolation basis on Gauss-Lobatto nodes.
%   This function calculates the values or derivatives of interpolation
%   basis functions defined on Gauss-Lobatto nodes within the interval
%   [-1/2, 1/2].
%
%   Syntax:
%      res = INTP_BASIS_GAUSS_LOBATTO(x, idx, deg)
%      res = INTP_BASIS_GAUSS_LOBATTO(x, idx, deg, drv)
%
%   Inputs:
%      x   - Scalar or array, evaluation points in [-1/2, 1/2].
%      idx - Index of the interpolation basis (1 <= idx <= deg+1).
%      deg - Degree of the interpolation basis (0 <= deg <= 3).
%      drv - (Optional) Derivative order (0 or 1). Default is 0.
%
%   Output:
%      res - Values of the interpolation basis or its derivative at `x`.
%
%   Example:
%      % Compute the degree-2 basis function for idx=1 at x=0
%      res = intp_basis_gauss_lobatto(0, 1, 2);
%
%      % Compute its derivative at x=0
%      res = intp_basis_gauss_lobatto(0, 1, 2, 1);
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.4

% Default argument
if nargin < 4
    drv = 0; % Default: no derivative
end

% Validate inputs
if deg < 0 || deg > 3
    error('Degree must be between 0 and 3.');
end
if idx < 1 || idx > deg + 1
    error('Index must be between 1 and the degree + 1.');
end
if drv < 0 || drv > 1
    error('Derivative order must be 0 or 1.');
end

% Dispatch to helper function
res = eval_basis(x, idx, deg, drv);
end

function res = eval_basis(x, idx, deg, drv)
% Dispatch based on degree
switch deg
    case 0
        res = eval_deg_0(x, idx, drv);
    case 1
        res = eval_deg_1(x, idx, drv);
    case 2
        res = eval_deg_2(x, idx, drv);
    case 3
        res = eval_deg_3(x, idx, drv);
end
end

% Degree 0
function res = eval_deg_0(x, idx, drv)
switch idx
    case 1
        res = eval_deg_0_1(x, drv);
end
end

function res = eval_deg_0_1(x, drv)
switch drv
    case 0
        res = ones(size(x));
    case 1
        res = zeros(size(x));
end
end

% Degree 1
function res = eval_deg_1(x, idx, drv)
switch idx
    case 1
        res = eval_deg_1_1(x, drv);
    case 2
        res = eval_deg_1_2(x, drv);
end
end

function res = eval_deg_1_1(x, drv)
switch drv
    case 0
        res = 1/2 - x;
    case 1
        res = -1 * ones(size(x));
end
end

function res = eval_deg_1_2(x, drv)
switch drv
    case 0
        res = 1/2 + x;
    case 1
        res = ones(size(x));
end
end

% Degree 2
function res = eval_deg_2(x, idx, drv)
switch idx
    case 1
        res = eval_deg_2_1(x, drv);
    case 2
        res = eval_deg_2_2(x, drv);
    case 3
        res = eval_deg_2_3(x, drv);
end
end

function res = eval_deg_2_1(x, drv)
switch drv
    case 0
        res = 2 * x .* (x - 1/2);
    case 1
        res = 4 * x - 1;
end
end

function res = eval_deg_2_2(x, drv)
switch drv
    case 0
        res = -4 * x .* x + 1;
    case 1
        res = -8 * x;
end
end

function res = eval_deg_2_3(x, drv)
switch drv
    case 0
        res = 2 * x .* (x + 1/2);
    case 1
        res = 4 * x + 1;
end
end

% Degree 3
function res = eval_deg_3(x, idx, drv)
switch idx
    case 1
        res = eval_deg_3_1(x, drv);
    case 2
        res = eval_deg_3_2(x, drv);
    case 3
        res = eval_deg_3_3(x, drv);
    case 4
        res = eval_deg_3_4(x, drv);
end
end

function res = eval_deg_3_1(x, drv)
switch drv
    case 0
        res = -x .* (x .* (5 * x - 5/2) - 1/4) - 1/8;
    case 1
        res = -15 * x .* x + 5 * x + 1/4;
end
end

function res = eval_deg_3_2(x, drv)
switch drv
    case 0
        res = -x .* (5 * sqrt(5)/4 - x .* (5 * sqrt(5) * x - 5/2)) + 5/8;
    case 1
        res = -15 * sqrt(5) * x .* x + 5 * x - 5 * sqrt(5)/4;
end
end

function res = eval_deg_3_3(x, drv)
switch drv
    case 0
        res = x .* (5 * sqrt(5)/4 - x .* (5 * sqrt(5) * x + 5/2)) + 5/8;
    case 1
        res = -15 * sqrt(5) * x .* x - 5 * x + 5 * sqrt(5)/4;
end
end

function res = eval_deg_3_4(x, drv)
switch drv
    case 0
        res = x .* (x .* (5 * x + 5/2) - 1/4) - 1/8;
    case 1
        res = 15 * x .* x + 5 * x - 1/4;
end
end
