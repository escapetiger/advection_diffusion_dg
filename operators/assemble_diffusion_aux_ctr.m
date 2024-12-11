function res = assemble_diffusion_aux_ctr(par, ref, dim, hx)
%ASSEMBLE_DIFFUSION_AUX_CTR Assemble diffusion matrix for an auxiliary 
% variable using central flux along a dimension.
%
%   This function computes the diffusion matrix for an auxiliary variable 
%   using central fluxes along the specified spatial dimension.
%
%   Syntax:
%      res = assemble_diffusion_aux_ctr(par, ref, dim, hx)
%
%   Inputs:
%      par - Parameters structure.
%      ref - Reference structure.
%      dim - Spatial dimension along which to compute the gradient.
%      hx  - Mesh size.
%
%   Outputs:
%      res - Sparse matrix.
%
%   Author: Yi Cai
%   Date: 2024-11-25
%   Version: 1.2

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = ref.n_dofs; % Number of local degrees of freedom per cell
ng = nl * nc; % Total number of global degrees of freedom
M = ref.v_u_vol; % Mass matrix
V = M \ ref.dv_u_vol{dim}; % Volume terms
FLi = M \ ref.dv_u_flx_i{2*dim-1}; % Left inflow flux terms
FLo = M \ ref.dv_u_flx_o{2*dim-1}; % Left outflow flux terms
FRi = M \ ref.dv_u_flx_i{2*dim}; % Right inflow flux terms
FRo = M \ ref.dv_u_flx_o{2*dim}; % Right outflow flux terms
h = 1 / hx(dim); % flux size

%========================================================================
% Initialize sparse matrix storage
%========================================================================
ns = nl^2; % Number of non-zero entries per cell
rows = zeros(ns*nc*5, 1);
cols = zeros(ns*nc*5, 1);
vals = zeros(ns*nc*5, 1);

%========================================================================
% Assemble sparse matrix
%========================================================================
% Volume integrals
a = 0;
[r, c, v] = find(-V*h);
i = 1:numel(v)*nc;
rows(a + i) = repmat(r, nc, 1) + kron((0:nc-1)' * nl, ones(numel(r), 1));
cols(a + i) = repmat(c, nc, 1) + kron((0:nc-1)' * nl, ones(numel(c), 1));
vals(a + i) = repmat(v, nc, 1);

% Right flux integrals
a = a + numel(i);
[r, c, v] = find(FRi*h/2);
i = 1:numel(v)*nc;
rows(a + i) = repmat(r, nc, 1) + kron((0:nc-1)' * nl, ones(numel(r), 1));
cols(a + i) = repmat(c, nc, 1) + kron((0:nc-1)' * nl, ones(numel(c), 1));
vals(a + i) = repmat(v, nc, 1);

a = a + numel(i);
[r, c, v] = find(FRo*h/2);
if par.bc(dim) == 0 % periodic B.C.
    m = multi_index(par.nx);
    m(:, dim) = m(:, dim) + 1;
    i = 1:numel(v)*nc;
    j = m2i(m, par.nx, par.bc(dim));
    rows(a + i) = repmat(r, nc, 1) + kron((0:nc-1)' * nl, ones(numel(r), 1));
    cols(a + i) = repmat(c, nc, 1) + kron((j-1) * nl, ones(numel(c), 1));
    vals(a + i) = repmat(v, nc, 1);
else % other B.C.
    m = multi_index(par.nx);
    m = m(m(:, dim) < par.nx(dim), :);
    ni = size(m, 1);
    i = 1:numel(v) * ni;
    j = m2i(m, par.nx, 1);
    rows(a+i) = repmat(r, ni, 1) + kron((j - 1)*nl, ones(numel(r), 1));
    m(:, dim) = m(:, dim) + 1;
    j = m2i(m, par.nx, 1);
    cols(a+i) = repmat(c, ni, 1) + kron((j - 1)*nl, ones(numel(c), 1));
    vals(a+i) = repmat(v, ni, 1);
end

% Left flux integrals
a = a + numel(i);
[r, c, v] = find(-FLi*h/2);
i = 1:numel(v)*nc;
rows(a + i) = repmat(r, nc, 1) + kron((0:nc-1)' * nl, ones(numel(r), 1));
cols(a + i) = repmat(c, nc, 1) + kron((0:nc-1)' * nl, ones(numel(c), 1));
vals(a + i) = repmat(v, nc, 1);

a = a + numel(i);
[r, c, v] = find(-FLo*h/2);
if par.bc(dim) == 0 % periodic B.C.
    m = multi_index(par.nx);
    m(:, dim) = m(:, dim) - 1;
    i = 1:numel(v) * nc;
    j = m2i(m, par.nx, 0);
    rows(a+i) = repmat(r, nc, 1) + kron((0:nc - 1)'*nl, ones(numel(r), 1));
    cols(a+i) = repmat(c, nc, 1) + kron((j - 1)*nl, ones(numel(c), 1));
    vals(a+i) = repmat(v, nc, 1);
else % other B.C.
    m = multi_index(par.nx);
    m = m(m(:, dim) > 1, :);
    ni = size(m, 1);
    i = 1:numel(v) * ni;
    j = m2i(m, par.nx, 1);
    rows(a+i) = repmat(r, ni, 1) + kron((j - 1)*nl, ones(numel(r), 1));
    m(:, dim) = m(:, dim) - 1;
    j = m2i(m, par.nx, 1);
    cols(a+i) = repmat(c, ni, 1) + kron((j - 1)*nl, ones(numel(c), 1));
    vals(a+i) = repmat(v, ni, 1);
end

%========================================================================
% Finalize sparse matrix
%========================================================================
nnz = find((rows > 0) & (cols > 0));
res = sparse(rows(nnz), cols(nnz), vals(nnz), ng, ng);

end
