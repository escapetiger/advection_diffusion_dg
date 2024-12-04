function [i] = m2i(m, s, bc, rc)
%M2I Convert a multi-dimensional index to a flat index (one-based).
%   This function converts a multi-dimensional index `m` into a one-based
%   flat linear index `i`, based on the array's shape `s` and its storage
%   order. It supports strict and circular boundary conditions.
%
%   Syntax:
%      i = M2I(m, s, bc)
%      i = M2I(m, s, bc, rc)
%
%   Inputs:
%      m  - Vector representing the multi-dimensional index.
%      s  - Vector specifying the maximal index in each dimension.
%      bc - Scalar specifying the boundary condition:
%           * 1: Strict boundary condition (default). Returns `i = 0`
%                for out-of-bound indices.
%           * 0: Circular boundary condition. Out-of-bound indices
%                wrap around using modulo arithmetic.
%      rc - (Optional) Storage order:
%           * 0 (default) for column-major order (MATLAB default).
%           * 1 for row-major order.
%
%   Outputs:
%      i  - Scalar representing the one-based flat linear index. Returns 0
%           if `m` is out of bounds with strict boundary conditions.
%
%   Example:
%      % Convert a valid index
%      m = [2, 3];
%      s = [4, 5];
%      i = m2i(m, s, 1, 0); % Column-major, strict boundary
%      disp(i); % Outputs 11
%
%      % Handle an out-of-bounds index with strict boundary
%      m = [5, 3];
%      i = m2i(m, s, 1, 0);
%      disp(i); % Outputs 0
%
%      % Handle an out-of-bounds index with circular boundary
%      m = [5, 3];
%      i = m2i(m, s, 0, 0);
%      disp(i); % Outputs a valid index
%
%   See also: IND2SUB, SUB2IND
%
%   Author: Yi Cai
%   Email: yicaim@stu.xmu.edu.cn
%   Date: 2024-11-21
%   Version: 1.1

if nargin < 4
    rc = 0; % Default to column-major storage order
end

if nargin < 3
    bc = 1; % Default to strict boundary condition
end

if size(m, ndims(m)) ~= numel(s)
    error('Multi-index and shape must have the same number of dimensions.');
end

if bc == 0
    m = mod(m - 1, s) + 1;
end

s = s';
p = ones(size(s));
if rc == 0
    p(2:end) = cumprod(s(1:end-1), 1);
else
    p(1:end-1) = flip(cumprod(flip(s(2:end), 1), 1), 1);
end
s = s';

i = (m - 1) * p + 1;

if bc == 1
    i(all(m > s, 2) | all(m < 1, 2)) = 0;
end

end
