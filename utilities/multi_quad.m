function [x, w] = multi_quad(n, a, b, quad_type, norm)
%MULTI_QUAD Generate multi-dimensional quadrature nodes and weights.
%   This function computes quadrature nodes and weights for a d-dimensional
%   domain using tensor-product Gaussian quadrature.
%
%   Syntax:
%      [x, w] = MQUAD(n, a, b, quad_type, norm)
%
%   Inputs:
%      n         - Number of quadrature nodes per dimension.
%      a         - Vector of lower bounds for each dimension.
%      b         - Vector of upper bounds for each dimension.
%      quad_type - Quadrature type:
%                  1 = Gauss-Legendre
%                  2 = Gauss-Lobatto
%      norm      - Logical, true to normalize weights. Default is false.
%
%   Outputs:
%      x - Matrix of quadrature nodes. Each column corresponds to one
%          dimension, and rows are the multi-dimensional nodes.
%      w - Row vector of corresponding quadrature weights.
%
%   Example:
%      % 2D quadrature with Gauss-Legendre
%      [x, w] = multi_quad(3, [-1, -1], [1, 1], 1, false);
%
%      % 3D quadrature with normalized Gauss-Lobatto
%      [x, w] = multi_quad(4, [0, 0, 0], [1, 1, 1], 2, true);
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

    if nargin < 5
        norm = false;
    end

    % Validate inputs
    if length(a) ~= length(b)
        error("Lower and upper bound vectors must have the same length.");
    end
    d = length(a); % Dimensionality of the domain

    if quad_type < 1 || quad_type > 2
        error("Invalid quadrature type. Use 1 (Gauss-Legendre) or 2 (Gauss-Lobatto).");
    end

    % Define quadrature functions
    quad_functions = {@gauss_legendre, @gauss_lobatto};

    % Generate nodes and weights for each dimension
    nodes = cell(1, d);
    weights = cell(1, d);
    for i = 1:d
        [nodes{i}, weights{i}] = quad_functions{quad_type}(n, a(i), b(i), norm);
    end

    % Create tensor product of nodes and weights
    [gridNodes{1:d}] = ndgrid(nodes{:});
    x = reshape(cat(d + 1, gridNodes{:}), [], d).'; % Multi-dimensional nodes

    [weightGrid{1:d}] = ndgrid(weights{:});
    w = prod(cat(d + 1, weightGrid{:}), d + 1); % Multi-dimensional weights
    w = w(:).'; % Row vector
end
