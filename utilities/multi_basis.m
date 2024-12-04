function [res] = multi_basis(varargin)
%MULTI_BASIS Calculate values of multi-dimensional polynomial basis.
%   This function computes the values or derivatives of a multi-dimensional
%   polynomial basis given the input points, degree, basis type, and
%   polynomial type.
%
%   Syntax:
%      res = MBASIS(x, deg, basis_t, poly_t)
%      res = MBASIS(x, deg, drv, drv_dim, basis_t, poly_t)
%
%   Inputs:
%      x        - Input points, where the first dimension corresponds to
%                 spatial dimensions.
%      deg      - Degree of the polynomial basis (non-negative integer).
%      drv      - (Optional) Order of the derivative (0 or 1). Default is 0.
%      drv_dim  - (Optional) Dimension for which to take the derivative.
%      basis_t  - Basis type:
%                 1 = Orthogonal Legendre basis
%                 2 = Interpolatory basis on Gauss-Legendre nodes
%                 3 = Interpolatory basis on Gauss-Lobatto nodes
%      poly_t   - Polynomial type:
%                 'P' for simplex polynomials
%                 'Q' for tensor-product polynomials
%
%   Outputs:
%      res - Array of basis function values. The first dimension
%            corresponds to the basis functions, while the remaining
%            dimensions correspond to the evaluation points.
%
%   Example:
%      % Compute values for Legendre polynomials
%      x = [-0.5, 0, 0.5];
%      res = multi_basis(x, 2, 1, 'P');
%
%      % Compute derivatives for Gauss-Lobatto basis
%      x = [-0.5, 0, 0.5];
%      res = multi_basis(x, 2, 1, 1, 3, 'Q');
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

% Handle inputs based on usage
if nargin == 4
    [x, deg, basis_t, poly_t] = deal(varargin{:});
    drv = 0; % Default derivative order
    drv_dim = []; % Irrelevant since drv = 0
elseif nargin == 6
    [x, deg, drv, drv_dim, basis_t, poly_t] = deal(varargin{:});
else
    error('Invalid number of inputs. Use either 4 or 6 arguments.');
end

% Validate inputs
sz = size(x);
n_dims = sz(1); % Number of spatial dimensions
eval_shape = sz(2:end); % Shape of evaluation points

if drv ~= 0 && drv ~= 1
    error('Derivative order "drv" must be 0 or 1.');
end
if drv == 1 && (drv_dim < 1 || drv_dim > n_dims)
    error('drv_dim must be a valid dimension index.');
end

% Generate multi-indices for polynomial type
if strcmp(poly_t, 'P')
    m = l1_index(n_dims, deg + 1); % Multi-indices for simplex polynomials
elseif strcmp(poly_t, 'Q')
    m = linf_index(n_dims, deg + 1); % Multi-indices for tensor-product polynomials
else
    error("Polynomial type must be 'P' or 'Q'.");
end

% Select appropriate basis function
basis_functions = {
    @orth_basis_legendre, ...
    @intp_basis_gauss_legendre, ...
    @intp_basis_gauss_lobatto
    };

if basis_t < 1 || basis_t > numel(basis_functions)
    error("Invalid basis type. Basis must be 1, 2, or 3.");
end
basis_fn = basis_functions{basis_t};

% Initialize output array
n_basis = size(m, 1); % Number of basis functions
res = ones([n_basis, eval_shape]);

% Compute basis values
for j = 1:n_basis
    for d = 1:n_dims
        idx_basis = repmat({':'}, 1, ndims(eval_shape)); idx_basis{1} = j;
        idx_dim = repmat({':'}, 1, ndims(eval_shape)); idx_dim{1} = d;

        % Determine if derivative is applied
        drv_mask = (drv == 1) && (d == drv_dim);

        % Evaluate basis function
        if basis_t == 1
            % Orthogonal Legendre basis
            res(idx_basis{:}) = res(idx_basis{:}) .* ...
                basis_fn(x(idx_dim{:}), m(j, d) - 1, drv_mask);
        else
            % Interpolatory basis
            res(idx_basis{:}) = res(idx_basis{:}) .* ...
                basis_fn(x(idx_dim{:}), m(j, d), deg, drv_mask);
        end
    end
end
end
