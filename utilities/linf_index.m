function [res] = linf_index(dim, n, rc)
%LINF_INDEX Generate one-based multi-indices with a uniform bound.
%   This function generates one-based multi-indices in `dim` dimensions
%   such that all indices satisfy the constraint:
%      max(index components) <= n.
%
%   Syntax:
%      res = LINF_INDEX(dim, n)
%      res = LINF_INDEX(dim, n, rc)
%
%   Inputs:
%      dim  - Number of dimensions (positive integer).
%      n    - Maximum value for any index component (positive integer).
%      rc   - (Optional) Storage order:
%             * 0 (default) for column-major (MATLAB default).
%             * 1 for row-major.
%
%   Output:
%      res  - An array of size (n^dim, dim), where each row is a multi-index
%             satisfying the constraint.
%
%   Example:
%      % Generate 2D indices with max index <= 3 in column-major order
%      res = linf_index(2, 3);
%
%      % Generate 3D indices with max index <= 2 in row-major order
%      res = linf_index(3, 2, 1);
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

    % Default storage order is column-major
    if nargin < 3
        rc = 0;
    end

    % Validate inputs
    if dim < 1 || n < 1 || floor(dim) ~= dim || floor(n) ~= n
        error('Both "dim" and "n" must be positive integers.');
    end

    % Special case for 1D
    if dim == 1
        res = (1:n).'; % Generate 1D indices directly
        return;
    end

    % Generate grid for multi-indices
    m = cell(1, dim);
    for d = 1:dim
        m{d} = 1:n; % Indices range from 1 to n
    end

    if rc == 0
        [m{1:dim}] = ndgrid(m{1:dim}); % Column-major
    else
        [m{1:dim}] = meshgrid(m{1:dim}); % Row-major
    end

    % Flatten the grids into rows of multi-indices
    res = zeros(n^dim, dim);
    for d = 1:dim
        res(:, d) = m{d}(:);
    end
end
