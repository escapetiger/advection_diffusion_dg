function res = assemble_diffusion_aux_bdr_ext_pos(par, dg, dim, hx)
%ASSEMBLE_DIFFUSION_AUX_BDR_EXT_POS Boundary effects of extraction 
% on the diffusion matrix for an auxiliary variable using positive flux.
%
%   This function computes the boundary part of sparse matrix for
%   an auxiliary variable introduced by extraction, when using positive
%   flux along the specific dimension.
%
%   Syntax:
%      res = assemble_diffusion_aux_bdr_ext_pos(par, dg, dim, hx)
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
FR = M \ dg.dv_u_flx_i{2*dim}; % Right inflow flux terms
h = 1 / hx(dim); % flux size

%========================================================================
% Assemble sparse matrix
%========================================================================
m = multi_index(par.nx);
m = m(m(:, dim) == par.nx(dim), :);
nb = size(m, 1);
i = m2i(m, par.nx, par.bc(dim) > 0);
[r, c, v] = find(FR*h);
rows = repmat(r, nb, 1) + kron((i-1) * nl, ones(numel(r), 1));
cols = repmat(c, nb, 1) + kron((i-1) * nl, ones(numel(c), 1));
vals = repmat(v, nb, 1);

%========================================================================
% Finalize sparse matrix
%========================================================================
nnz = find((rows > 0) & (cols > 0));
res = sparse(rows(nnz), cols(nnz), vals(nnz), ng, ng);

end
