function [res] = multi_index(s, rc)
%MULTI_INDEX Generate one-based multi-indices for a given shape.
%   This function generates one-based multi-indices for a given shape `s`.
%   Each row of the output contains a valid multi-index, satisfying:
%      1 <= index(dim) <= s(dim) for all dimensions.
%
%   Syntax:
%      res = MULTI_INDEX(s)
%      res = MULTI_INDEX(s, rc)
%
%   Inputs:
%      s  - Shape of the indices (positive integer vector).
%      rc - (Optional) Storage order:
%           * 0 (default) for column-major (MATLAB default).
%           * 1 for row-major.
%
%   Output:
%      res - An array of size (prod(s), numel(s)), where each row is a 
%            multi-index satisfying the shape constraint.
%
%   Example:
%      % Generate 2D indices for a shape [3, 2] in column-major order
%      res = multi_index([3, 2]);
%
%      % Generate 3D indices for a shape [2, 2, 2] in row-major order
%      res = multi_index([2, 2, 2], 1);
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.0

    % Default storage order is column-major
    if nargin < 2
        rc = 0;
    end

    % Validate inputs
    if ~isvector(s) || any(s < 1) || any(floor(s) ~= s)
        error('"s" must be a vector of positive integers.');
    end
    if ~isscalar(rc) || ~ismember(rc, [0, 1])
        error('"rc" must be 0 (column-major) or 1 (row-major).');
    end

    dim = numel(s);

    % Special case for 1D
    if dim == 1
        res = (1:s).'; % Generate 1D indices directly
        return;
    end

    % Generate grid for multi-indices
    m = cell(1, dim);
    for d = 1:dim
        m{d} = 1:s(d); % Indices range from 1 to s(d)
    end

    % Determine storage order
    if rc == 0
        [m{1:dim}] = ndgrid(m{1:dim}); % Column-major
    else
        [m{1:dim}] = meshgrid(m{1:dim}); % Row-major
    end

    % Flatten the grids into rows of multi-indices
    res = zeros(prod(s), dim);
    for d = 1:dim
        res(:, d) = m{d}(:);
    end
end
