function [res] = intp_basis_gauss_legendre(x, idx, deg, drv)
%INTP_BASIS_GAUSS_LEGENDRE Compute interpolation basis on Gauss-Legendre nodes.
%   This function calculates the values or derivatives of interpolation
%   basis functions defined on Gauss-Legendre nodes within [-1/2, 1/2].
%
%   Syntax:
%      res = INTP_BASIS_GAUSS_LEGENDRE(x, idx, deg)
%      res = INTP_BASIS_GAUSS_LEGENDRE(x, idx, deg, drv)
%
%   Inputs:
%      x   - Scalar or array, evaluation points in [-1/2, 1/2].
%      idx - Index of the interpolation basis (zero-based, 0 <= idx <= deg).
%      deg - Degree of the interpolation basis (0 <= deg <= 3).
%      drv - (Optional) Derivative order (0 or 1). Default is 0.
%
%   Output:
%      res - Values of the interpolation basis or its derivative at `x`.
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

if nargin < 4
    drv = 0; % Default derivative order is 0
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

% Dispatch to degree-specific helper function
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
        res = eval_0_1(x, drv);
end
end

function res = eval_0_1(x, drv)
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
        res = -sqrt(3) * x + 1 / 2;
    case 1
        res = -sqrt(3) * ones(size(x));
end
end

function res = eval_deg_1_2(x, drv)
switch drv
    case 0
        res = sqrt(3) * x + 1 / 2;
    case 1
        res = sqrt(3) * ones(size(x));
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
        res = x .* (10 / 3 * x - sqrt(5) / sqrt(3));
    case 1
        res = 20 / 3 * x - sqrt(5) / sqrt(3);
end
end

function res = eval_deg_2_2(x, drv)
switch drv
    case 0
        res = 1 - 20 / 3 * x.^2;
    case 1
        res = -40 / 3 * x;
end
end

function res = eval_deg_2_3(x, drv)
switch drv
    case 0
        res = x .* (10 / 3 * x + sqrt(5) / sqrt(3));
    case 1
        res = 20 / 3 * x + sqrt(5) / sqrt(3);
end
end

% Degree 3
function res = eval_deg_3(x, idx, drv)
a = 0.430568155797026;
b = 0.169990521792428;
c = (a + b) * (a - b);
switch idx
    case 1
        res = eval_deg_3_1(x, drv, a, b, c);
    case 2
        res = eval_deg_3_2(x, drv, a, b, c);
    case 3
        res = eval_deg_3_3(x, drv, a, b, c);
    case 4
        res = eval_deg_3_4(x, drv, a, b, c);
end
end

function res = eval_deg_3_1(x, drv, a, b, c)
switch drv
    case 0
        res = x .* (x .* (1 / (2 * c) - x / (2 * a * c)) - ...
            b * (1 / c - b / (a * c)) / 2 + b / (2 * c)) - b^2 / (2 * c);
    case 1
        res = x .* (1 / c - (3 * x) / (2 * a * c)) - ...
            (b * (1 / c - b / (a * c))) / 2 + b / (2 * c);
end
end

function res = eval_deg_3_2(x, drv, a, b, c)
switch drv
    case 0
        res = a^2 / (2 * c) - x .* (x .* (1 / (2 * c) - x / (2 * b * c)) + ...
            a^2 / (2 * b * c));
    case 1
        res = -x .* (1 / c - (3 * x) / (2 * b * c)) - a^2 / (2 * b * c);
end
end

function res = eval_deg_3_3(x, drv, a, b, c)
switch drv
    case 0
        res = a^2 / (2 * c) - x .* (x .* ((1 / (a + b) + a / (b * (a + b))) / ...
            (2 * (a - b)) - a / (2 * b * c) + x / (2 * b * c)) + a / (2 * c) ...
            - a * (1 / (a + b) + a / (b * (a + b))) / (2 * (a - b)));
    case 1
        res = (a * (1 / (a + b) + a / (b * (a + b)))) / (2 * (a - b)) - a / ...
            (2 * c) - x .* ((1 / (a + b) + a / (b * (a + b))) / (a - b) - ...
            a / (b * c) + (3 * x) / (2 * b * c));
end
end

function res = eval_deg_3_4(x, drv, a, b, c)
switch drv
    case 0
        res = x .* (x .* ((1 / (a + b) + b / (a * (a + b))) / (2 * (a - b)) - ...
            b / (2 * a * c) + x / (2 * a * c)) + b / (2 * c) - ...
            (b * (1 / (a + b) + b / (a * (a + b)))) / (2 * (a - b))) - ...
            b^2 / (2 * c);
    case 1
        res = x .* ((1 / (a + b) + b / (a * (a + b))) / (a - b) - b / (a * c) + ...
            (3 * x) / (2 * a * c)) + b / (2 * c) - ...
            (b * (1 / (a + b) + b / (a * (a + b)))) / (2 * (a - b));
end
end
