function res = assemble_diffusion_aux_bdr_jmp_ctr(par, dg, dim, hx)
%ASSEMBLE_DIFFUSION_AUX_BDR_JMP_CTR Boundary effects of jump 
% on the diffusion matrix for an auxiliary variable using central flux.
% 
%   Syntax:
%      res = assemble_diffusion_aux_bdr_jmp_ctr(par, dg, dim, hx)
%
%   Inputs:
%      par - Parameters structure containing.
%      dg  - DG structure containing.
%      dim - Spatial dimension along which to compute the gradient.
%      hx  - Mesh size.
%
%   Outputs:
%      res - Sparse matrix.
%
%   Author: Yi Cai
%   Date: 2024-11-28
%   Version: 1.1

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = dg.n_dofs; % Number of local degrees of freedom per cell
ng = nl * nc; % Total number of global degrees of freedom
M = dg.v_u_vol; % Mass matrix
FL = M \ dg.dv_u_flx_i{2*dim-1}; % Left inflow flux terms
FR = M \ dg.dv_u_flx_i{2*dim}; % Right inflow flux terms
h = 1 / hx(dim); % flux size
ns = nl^2; % Number of non-zero entries per cell

%========================================================================
% Assemble sparse matrix
%========================================================================
m = multi_index(par.nx);
m1 = m(m(:, dim) == 1, :);
m2 = m(m(:, dim) == par.nx(dim), :);
nb1 = size(m1, 1);
nb2 = size(m2, 1);
rows = zeros(ns*(nb1+nb2), 1);
cols = zeros(ns*(nb1+nb2), 1);
vals = zeros(ns*(nb1+nb2), 1);

a = 0;
[r, c, v] = find(-FL*h/2);
i = m2i(m1, par.nx, par.bc(dim) > 0);
j = 1:numel(v)*nb1;
rows(a + j) = repmat(r, nb1, 1) + kron((i-1) * nl, ones(numel(r), 1));
cols(a + j) = repmat(c, nb1, 1) + kron((i-1) * nl, ones(numel(c), 1));
vals(a + j) = repmat(v, nb1, 1);

a = a + numel(j);
[r, c, v] = find(-FR*h/2);
i = m2i(m2, par.nx, par.bc(dim) > 0);
j = 1:numel(v)*nb2;
rows(a + j) = repmat(r, nb2, 1) + kron((i-1) * nl, ones(numel(r), 1));
cols(a + j) = repmat(c, nb2, 1) + kron((i-1) * nl, ones(numel(c), 1));
vals(a + j) = repmat(v, nb2, 1);

%========================================================================
% Finalize sparse matrix
%========================================================================
nnz = find((rows > 0) & (cols > 0));
res = sparse(rows(nnz), cols(nnz), vals(nnz), ng, ng);

end
