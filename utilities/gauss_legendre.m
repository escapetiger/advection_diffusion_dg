function [x, w] = gauss_legendre(n, a, b, norm)
%GAUSS_LEGENDRE Generate Gauss-Legendre quadrature nodes and weights.
%   This function computes the nodes (x) and weights (w) for the
%   Gauss-Legendre quadrature rule, which integrates polynomials of degree
%   up to 2n-1 exactly. The nodes and weights can be computed within the
%   interval [a, b] and optionally normalized.
%
%   Syntax:
%      [x, w] = GAUSS_LEGENDRE(n)
%      [x, w] = GAUSS_LEGENDRE(n, a, b)
%      [x, w] = GAUSS_LEGENDRE(n, a, b, norm)
%
%   Inputs:
%      n    - Number of quadrature nodes (positive integer, n >= 1).
%      a    - (Optional) Lower bound of the interval [a, b]. Default is -1.
%      b    - (Optional) Upper bound of the interval [a, b]. Default is 1.
%      norm - (Optional) Logical flag for normalizing weights:
%             * false (default): Use standard quadrature weights.
%             * true: Normalize weights such that they sum to 1.
%
%   Outputs:
%      x - Row vector of n quadrature nodes within [a, b].
%      w - Row vector of n quadrature weights corresponding to the nodes.
%
%   Example:
%      % Compute 5 Gauss-Legendre nodes and weights on [-1, 1]
%      [x, w] = gauss_legendre(5);
%      disp(x); % Quadrature nodes
%      disp(w); % Quadrature weights
%
%      % Compute 5 Gauss-Legendre nodes and weights on [0, 2]
%      [x, w] = gauss_legendre(5, 0, 2, false);
%
%   See also: LEGENDRE, EIG
%
%   Author: Yi Cai
%   Email: yicaim@stu.xmu.edu.cn
%   Date: 2024-11-21
%   Version: 1.1

% Assign default values if not provided
if nargin < 4
    norm = false; % Default: Not normalize weights
end
if nargin < 3
    a = -1; b = 1; % Default interval: [-1, 1]
end

% Special case: a == b, integral reduces to a single point
if a == b
    x = a;
    w = 1;
    return;
end

% Compute Gauss-Legendre nodes and weights on [-1, 1]
beta = 0.5 * (1 - (2 * (1:n-1)).^(-2)).^(-0.5); % Recurrence coefficients
T = diag(beta, 1) + diag(beta, -1); % Jacobi matrix
[V, D] = eig(T); % Eigenvalue decomposition
x = diag(D)'; % Nodes (eigenvalues)
w = 2 * (V(1, :).^2); % Weights

% Map nodes and weights from [-1, 1] to [a, b]
x = (a * (1 - x) + b * (1 + x)) / 2; % Linear map
w = (b - a) / 2 * w; % Scale weights

% Normalize weights if requested
if norm
    w = w / (b - a);
end

end
