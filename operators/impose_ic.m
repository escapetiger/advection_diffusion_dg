function res = impose_ic(par, ref, x, hx)
%IMPOSE_IC Impose initial condition for DG simulation.
%
%   This function computes the initial condition for a DG method using
%   projection onto the basis functions.
%
%   Syntax:
%      res = impose_ic(par, ref, x, hx)
%
%   Inputs:
%      par - Parameters structure.
%      ref  - DG structure.
%      x   - Cell array of grid points in each dimension.
%      hx  - Grid spacing in each dimension.
%
%   Outputs:
%      res - Initial condition projected onto the DG basis.
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = ref.n_dofs; % Number of local degrees of freedom
ng = nl * nc; % Total number of global degrees of freedom
nq = ref.nq_vol; % Number of quadrature points per cell
xq = ref.xq_vol; % Quadrature points (reference)
wq = ref.wq_vol; % Quadrature weights (reference)
vq = ref.vq_vol; % Basis function values at quadrature points
M = ref.v_u_vol; % Mass matrix for volume

%========================================================================
% Points to evaluate
%========================================================================
m = multi_index(par.nx);
X = cell(1, par.dim);
for d = 1:par.dim
    X{d} = reshape(x{d}(m(:, d)), 1, []) + xq(d, :)' * hx(d);
end

%========================================================================
% Projection/Interpolation
%========================================================================
if par.basis_t == 1 || par.poly_t == 'P'
    U = reshape(par.fn_ic(par, X), nq, nc);
    P = reshape(vq, nl, nq);
    W = diag(wq);
    res = reshape(M \ (P * W * U), ng, 1);
else
    res = par.fn_ic(par, X);
    res = res(:);
end

end
