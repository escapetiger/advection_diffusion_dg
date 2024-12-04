function [x, w] = gauss_lobatto(n, a, b, norm)
%GAUSS_LOBATTO Generate Gauss-Lobatto quadrature nodes and weights.
%   This function computes the nodes (x) and weights (w) for the 
%   Gauss-Lobatto quadrature rule. The rule includes the endpoints of 
%   the interval and is useful for numerical integration.
%
%   Syntax:
%      [x, w] = GAUSS_LOBATTO(n)
%      [x, w] = GAUSS_LOBATTO(n, a, b)
%      [x, w] = GAUSS_LOBATTO(n, a, b, norm)
%
%   Inputs:
%      n    - Number of quadrature nodes (n >= 2).
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
%      % Compute 5 Gauss-Lobatto nodes and weights on [-1, 1]
%      [x, w] = gauss_lobatto(5);
%      disp(x); % Quadrature nodes
%      disp(w); % Quadrature weights
%
%      % Compute 5 Gauss-Lobatto nodes and weights on [0, 2]
%      [x, w] = gauss_lobatto(5, 0, 2, false);
%
%   See also: LEGENDRE, GAUSS_LEGENDRE
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

    % Initialize nodes with Chebyshev-Gauss-Lobatto points
    x = cos(pi * ((n-1):-1:0) / (n-1))'; % Initial guess
    P = zeros(n, n); % Legendre polynomial matrix
    xold = 2 * ones(size(x)); % Initial placeholder for convergence check

    % Iteratively refine nodes using Newton's method
    while max(abs(x - xold)) > eps(1)
        xold = x;

        % Compute Legendre polynomials at x
        P(1, :) = 1; % P_0(x)
        P(2, :) = x'; % P_1(x)
        for k = 3:n
            P(k, :) = ((2*k-3)*x'.*P(k-1, :) - (k-2)*P(k-2, :)) / (k-1);
        end

        % Newton's method update
        x = xold - (x .* P(n, :)' - P(n-1, :)') ./ (n * P(n, :)');
    end

    % Compute weights
    w = 2 ./ ((n-1)*n * P(n, :)'.^2);

    % Map nodes and weights from [-1, 1] to [a, b]
    x = (a * (1 - x) + b * (1 + x)) / 2; % Linear map
    w = (b - a) / 2 * w; % Scale weights
    x = x'; w = w';

    % Normalize weights if requested
    if norm
        w = w / (b - a);
    end
end
