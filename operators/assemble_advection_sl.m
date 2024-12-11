function res = assemble_advection_sl(par, ref, kt, kv)
%ASSEMBLE_ADVECTION_SL Assemble advection matrix using SLDG approach.
%
%   This function computes the advection matrix using the SL method.
%
%   Syntax:
%      res = assemble_advection_sl(par, ref, kt, kv)
%
%   Inputs:
%      par  - Parameters structure.
%      ref  - Reference structure.
%      hx   - Mesh size.
%      kt   - Time index.
%      kv   - Velocity index.
%
%   Outputs:
%      res - Sparse matrix.
%
%   Author: Yi Cai
%   Date: 2024-12-8
%   Version: 1.1

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = ref.n_dofs; % Number of local degrees of freedom per cell
ng = nl * nc; % Total number of global degrees of freedom
ktv = m2i([kt, kv], [ref.nt, ref.nv]); % Linear index of samples
np = ref.n_volps(ktv); % Total number of volume pieces
M = arrayfun( ...
    @(k) ref.v_u_vol \ ref.v_u_volp{ktv}(:, :, k), 1:np, ...
    'Un', false); % Mass piece
s = ref.volp_shift{ktv}; % Index shifts

%========================================================================
% Assembling
%========================================================================
ns = nl^2; % Number of non-zero entries per cell
rows = zeros(ns*nc*np, 1);
cols = zeros(ns*nc*np, 1);
vals = zeros(ns*nc*np, 1);

a = 0;
for p = 1:np
    [r, c, v] = find(M{p});
    if all(par.bc == 0) % periodic B.C.
        m = multi_index(par.nx);
        i = 1:numel(v) * nc;
        j = m2i(m + s(:, p)', par.nx, 0) - 1;
        k = (0:nc-1)';
        rows(a+i) = repmat(r, nc, 1) + kron(k*nl, ones(numel(r), 1));
        cols(a+i) = repmat(c, nc, 1) + kron(j*nl, ones(numel(c), 1));
        vals(a+i) = repmat(v, nc, 1);
    else
        m = multi_index(par.nx);
        bi = (m + s(:, p)' < 1) .* (2 * (1:par.dim) - 1) + ...
             (m + s(:, p)' > par.nx) .* (2 * (1:par.dim));
        m = m(sum(bi~=0, 2) == 0, :);
        ni = size(m, 1);
        i = 1:numel(v) * ni;
        j = m2i(m + s(:, p)', par.nx, 1) - 1;
        k = m2i(m, par.nx, 1) - 1;
        rows(a+i) = repmat(r, ni, 1) + kron(k*nl, ones(numel(r), 1));
        cols(a+i) = repmat(c, ni, 1) + kron(j*nl, ones(numel(c), 1));
        vals(a+i) = repmat(v, ni, 1);
    end
    a = a + numel(i);
end

%========================================================================
% Finalize sparse matrix
%========================================================================
nnz = find((rows > 0) & (cols > 0));
res = sparse(rows(nnz), cols(nnz), vals(nnz), ng, ng);

end
