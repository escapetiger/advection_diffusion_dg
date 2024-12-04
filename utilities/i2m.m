function [m] = i2m(i, s, rc)
%I2M Convert a flat index to a multi-index (one-based).
%   This function converts a one-based flat linear index `i` into a 
%   multi-dimensional index `m` based on the array's shape `s` and its 
%   storage order.
%
%   Syntax:
%      m = I2M(i, s)
%      m = I2M(i, s, rc)
%
%   Inputs:
%      i  - Vector of one-based flat linear indices.
%      s  - Vector specifying the shape of the array (dimensions).
%      rc - (Optional) Storage order:
%           * 0 (default) for column-major order (MATLAB default).
%           * 1 for row-major order.
%
%   Outputs:
%      m  - Matrix where each row corresponds to the multi-dimensional
%           index of the corresponding element in `i`.
%
%   Example:
%      % Convert a valid flat index to a multi-dimensional index
%      i = 11;
%      s = [4, 5];
%      m = i2m(i, s, 0); % Column-major
%      disp(m); % Outputs [3, 3]
%
%   See also: M2I
%
%   Author: Yi Cai
%   Email: yicaim@stu.xmu.edu.cn
%   Date: 2024-11-21
%   Version: 1.1

    if nargin < 3
        rc = 0; % Default to column-major order
    end

    % Validate input
    if any((i > prod(s)) | (i < 1))
        error('Flat index out of range.');
    end

    % Reshape s to ensure it is a column vector
    s = reshape(s, [], 1);

    % Compute strides for storage order
    p = ones(numel(s), 1);
    if rc == 0 % Column-major order
        p(2:end) = cumprod(s(1:end-1), 1);
    else % Row-major order
        p(1:end-1) = flip(cumprod(flip(s(2:end), 1), 1), 1);
    end

    % Reshape inputs for consistent operations
    i = reshape(i, [], 1); % Ensure i is a column vector
    p = reshape(p, 1, []); % Ensure p is a row vector
    s = reshape(s, 1, []); % Ensure s is a row vector

    % Compute the multi-dimensional index
    m = mod(floor((i - 1) ./ p), s) + 1;
end
