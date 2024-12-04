function [res] = l1_index(dim, n)
%L1_INDEX Generate one-based multi-indices with limited l1-norm.
%   This function generates one-based multi-indices in `dim` dimensions
%   such that their l1-norm is no more than `n`.
%
%   Syntax:
%      res = L1_INDEX(dim, n)
%
%   Inputs:
%      dim  - Number of dimensions (positive integer).
%      n    - Maximum l1-norm of the indices (non-negative integer).
%
%   Output:
%      res  - An array of size (N, dim), where N is the number of
%             multi-indices satisfying the constraint.
%
%   Example:
%      % Generate 2D indices with l1-norm <= 2
%      res = l1_index(2, 2);
%
%      % Generate 3D indices with l1-norm <= 1
%      res = l1_index(3, 1);
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.2

    % Validate inputs
    if dim < 1 || n < 0 || floor(dim) ~= dim || floor(n) ~= n
        error('Both "dim" and "n" must be non-negative integers, with dim >= 1.');
    end

    % Special case for 1D
    if dim == 1
        res = (1:n).'; % One-based indices
        return;
    end

    % Generate multi-indices for all l1-norms up to `n`
    res = cell(n, 1);
    for k = 0:n-1
        % Generate dividers for partitioning `k` into `dim` parts
        dividers = flipud(nchoosek(1:(k + dim - 1), dim - 1));
        res{k + 1} = [dividers(:, 1), diff(dividers, 1, 2), ...
            (k + dim) - dividers(:, end)];
    end

    % Concatenate results into a single array
    res = cell2mat(res);
end
