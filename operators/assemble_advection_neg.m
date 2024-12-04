function res = assemble_advection_neg(par, dg, dim, hx)
%ASSEMBLE_ADVECTION_NEG Assemble advection matrix using negative flux 
% along a dimension.
%
%   This function computes the advection matrix for negative fluxes along
%   the specified spatial dimension.
%
%   Syntax:
%      res = assemble_advection_neg(par, dg, dim, hx)
%
%   Inputs:
%      par - Parameters structure.
%      dg  - DG structure.
%      dim - Spatial dimension along which to compute the gradient.
%      hx  - Mesh size.
%
%   Outputs:
%      res - Sparse gradient matrix for negative flux along the specified dimension.
%
%   Author: Yi Cai
%   Date: 2024-11-25
%   Version: 1.2

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = dg.n_dofs; % Number of local degrees of freedom per cell
ng = nl * nc; % Total number of global degrees of freedom
M = dg.v_u_vol; % Mass matrix
V = M \ dg.dv_u_vol{dim}; % Volume terms
FL = M \ dg.dv_u_flx_o{2*dim-1}; % Left flux terms
FR = M \ dg.dv_u_flx_i{2*dim}; % Right flux terms
h = 1 / hx(dim); % flux size

%========================================================================
% Initialize sparse matrix storage
%========================================================================
ns = nl^2; % Number of non-zero entries per cell
rows = zeros(ns*nc*3, 1);
cols = zeros(ns*nc*3, 1);
vals = zeros(ns*nc*3, 1);

%========================================================================
% Assemble sparse matrix
%========================================================================
% Volume integrals
a = 0;
[r, c, v] = find(-V*h);
i = 1:numel(v) * nc;
rows(a+i) = repmat(r, nc, 1) + kron((0:nc - 1)'*nl, ones(numel(r), 1));
cols(a+i) = repmat(c, nc, 1) + kron((0:nc - 1)'*nl, ones(numel(c), 1));
vals(a+i) = repmat(v, nc, 1);

% Right flux integrals
a = a + numel(i);
[r, c, v] = find(FR*h);
i = 1:numel(v) * nc;
rows(a+i) = repmat(r, nc, 1) + kron((0:nc - 1)'*nl, ones(numel(r), 1));
cols(a+i) = repmat(c, nc, 1) + kron((0:nc - 1)'*nl, ones(numel(c), 1));
vals(a+i) = repmat(v, nc, 1);

% Left flux integrals
a = a + numel(i);
[r, c, v] = find(-FL*h);
m = multi_index(par.nx);
m(:, dim) = m(:, dim) - 1;
i = 1:numel(v) * nc;
j = m2i(m, par.nx, 0);
rows(a+i) = repmat(r, nc, 1) + kron((0:nc - 1)'*nl, ones(numel(r), 1));
cols(a+i) = repmat(c, nc, 1) + kron((j - 1)*nl, ones(numel(c), 1));
vals(a+i) = repmat(v, nc, 1);

%========================================================================
% Finalize sparse matrix
%========================================================================
nnz = find((rows > 0) & (cols > 0));
res = sparse(rows(nnz), cols(nnz), vals(nnz), ng, ng);

end
