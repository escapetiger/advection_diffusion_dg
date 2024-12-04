function [x, w] = gauss_chebyshev(n, a, b, norm)
%GAUSS_CHEBYSHEV Generate Gauss-Chebyshev quadrature nodes and weights.
%   This function computes the nodes (x) and weights (w) for the 
%   Gauss-Chebyshev quadrature rule, which is exact for polynomials 
%   of degree up to 2n-1 under the weight function (1-x^2)^(-1/2).
%
%   Syntax:
%      [x, w] = GAUSS_CHEBYSHEV(n)
%      [x, w] = GAUSS_CHEBYSHEV(n, a, b)
%      [x, w] = GAUSS_CHEBYSHEV(n, a, b, norm)
%
%   Inputs:
%      n    - Number of quadrature nodes (n >= 1).
%      a    - (Optional) Lower bound of the interval [a, b]. Default is -1.
%      b    - (Optional) Upper bound of the interval [a, b]. Default is 1.
%      norm - (Optional) Logical flag for normalizing weights:
%             * false (default): Use standard Gauss-Chebyshev weights.
%             * true: Normalize weights such that they sum to 1.
%
%   Outputs:
%      x - Row vector of n quadrature nodes within [a, b].
%      w - Row vector of n quadrature weights corresponding to the nodes.
%
%   Example:
%      % Compute 5 Gauss-Chebyshev nodes and weights on [-1, 1]
%      [x, w] = gauss_chebyshev(5);
%      disp(x); % Quadrature nodes
%      disp(w); % Quadrature weights
%
%      % Compute 5 Gauss-Chebyshev nodes and weights on [0, 2]
%      [x, w] = gauss_chebyshev(5, 0, 2, false);
%
%   See also: LEGENDRE, GAUSS_LOBATTO
%
%   Author: Yi Cai
%   Email: yicaim@stu.xmu.edu.cn
%   Date: 2024-11-21
%   Version: 1.1

    % Assign default values if not provided
    if nargin < 4
        norm = false; % Default: No normalization
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

    % Compute Chebyshev nodes and weights
    i = (1:n)';
    x = cos(pi * (i - 0.5) / n)'; % Chebyshev nodes on [-1, 1]
    w = pi / n * ones(1, n); % Uniform weights

    % Map nodes and weights from [0, 1] to [a, b]
    x = a + (b - a) * x;
    w = (b - a) * w;

    % Normalize weights if requested
    if norm
        w = w / (b - a);
    end
end
