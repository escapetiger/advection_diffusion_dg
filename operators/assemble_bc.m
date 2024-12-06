function res = assemble_bc(par, dg, x, hx, t, ht)
%ASSEMBLE_BC
%   Assemble vectors for boundary conditions of diffusion terms.

%========================================================================
% Parse parameters and initialize variables
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = dg.n_dofs; % Number of local DOFs per cell
ng = nl * nc; % Total number of global DOFs
nf = dg.n_flxs; % Number of fluxes
nq = dg.nq_flx; % Number of flux quadrature points
xq = dg.xq_flx; % Flux quadrature points
M = dg.v_u_vol; % Mass Matrix

res = cell(1, nf);
for i = 1:par.dim
   res{i} = sparse(ng, 1);
end

%========================================================================
% Apply boundary conditions for each dimension
%========================================================================
for d = 1:par.dim
    if par.bc(d) == 0, continue; end % omit periodic B.C.

    %--------------------------------------------------------------------
    % Left boundary
    %--------------------------------------------------------------------
    j = repmat({':'}, 1, par.dim);
    if par.dim == 1
        j{d} = 1;
    else
        j{d} = ones(length(x{d}), 1);
    end
    m = multi_index(par.nx);
    b = m(m(:, d) == 1, :);
    nb = size(b, 1);
    k = m2i(b, par.nx);
    rows = repmat((k - 1)' * nl, nl, 1) + repmat((1:nl)', 1, nb);
    rows = rows(:);
    cols = ones(size(rows));

    % Points to evaluate
    X = cell(1, par.dim);
    for dd = 1:par.dim
        X{dd} = reshape(x{dd}(b(:, dd)), 1, []) + xq{2*d-1}(dd, :)' * hx(dd);
        X{dd} = X{dd}(:);
    end

    % Evaluation
    F = par.fn_bc(par, X, t + ht, 2*d-1);
    F = reshape(F, nq, nb);
    P = dg.vq_flx{2*d-1};
    W = diag(dg.wq_flx{2*d-1});
    vals = M \ (P * W * F) / hx(d);
    vals = vals(:);
    res{2*d-1} = sparse(rows, cols, vals, ng, 1);

    %--------------------------------------------------------------------
    % Right boundary
    %--------------------------------------------------------------------
    j = repmat({':'}, 1, par.dim);
    if par.dim == 1
        j{d} = par.nx(d);
    else
        j{d} = ones(length(x{d}), 1) * par.nx(d);
    end
    m = multi_index(par.nx);
    b = m(m(:, d) == par.nx(d), :);
    nb = size(b, 1);
    k = m2i(b, par.nx, 1);
    rows = repmat((k - 1)' * nl, nl, 1) + repmat((1:nl)', 1, nb);
    rows = rows(:);
    cols = ones(size(rows));

    % Points to evaluate
    X = cell(1, par.dim);
    for dd = 1:par.dim
        X{dd} = reshape(x{dd}(b(:, dd)), 1, []) + xq{2*d}(dd, :)' * hx(dd);
    end

    % Evaluation
    F = par.fn_bc(par, X, t + ht, 2*d);
    F = reshape(F, nq, nb);
    P = dg.vq_flx{2*d};
    W = diag(dg.wq_flx{2*d});
    vals = M \ (P * W * F) / hx(d);
    vals = vals(:);
    res{2*d} = sparse(rows, cols, vals, ng, 1);
end

end
