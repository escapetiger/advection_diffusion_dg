function res = assemble_bc(par, ref, x, hx, t, ht)
%ASSEMBLE_BC
%   Assemble vectors for boundary conditions.

%========================================================================
% Parse parameters and initialize variables
%========================================================================
if all(par.bc == 0) % Nothing to do with periodic B.C.
    res = []; 
    return; 
end

nc = prod(par.nx); % Total number of cells
nl = ref.n_dofs; % Number of local DOFs per cell
ng = nl * nc; % Total number of global DOFs
nf = ref.n_flxs; % Number of fluxes
nq = ref.nq_flx; % Number of flux quadrature points
xq = ref.xq_flx; % Flux quadrature points
M = ref.v_u_vol; % Mass Matrix

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
    m = multi_index(par.nx);
    b = m(m(:, d) == 1, :);
    nb = size(b, 1);
    k = m2i(b, par.nx);
    rows = repmat((k - 1)'*nl, nl, 1) + repmat((1:nl)', 1, nb);
    rows = rows(:);
    cols = ones(size(rows));

    % Points to evaluate
    X = cell(1, par.dim);
    for dd = 1:par.dim
        X{dd} = reshape(x{dd}(b(:, dd)), 1, []) + xq{2*d-1}(dd, :)' * hx(dd);
        X{dd} = X{dd}(:);
    end

    % Evaluation
    F = par.fn_bc(par, X, t+ht, 2*d-1);
    F = reshape(F, nq, nb);
    P = ref.vq_flx{2*d-1};
    W = diag(ref.wq_flx{2*d-1});
    vals = M \ (P * W * F) / hx(d);
    vals = vals(:);
    res{2*d-1} = sparse(rows, cols, vals, ng, 1);

    %--------------------------------------------------------------------
    % Right boundary
    %--------------------------------------------------------------------
    m = multi_index(par.nx);
    b = m(m(:, d) == par.nx(d), :);
    nb = size(b, 1);
    k = m2i(b, par.nx, 1);
    rows = repmat((k - 1)'*nl, nl, 1) + repmat((1:nl)', 1, nb);
    rows = rows(:);
    cols = ones(size(rows));

    % Points to evaluate
    X = cell(1, par.dim);
    for dd = 1:par.dim
        X{dd} = reshape(x{dd}(b(:, dd)), 1, []) + xq{2*d}(dd, :)' * hx(dd);
        X{dd} = X{dd}(:);
    end

    % Evaluation
    F = par.fn_bc(par, X, t+ht, 2*d);
    F = reshape(F, nq, nb);
    P = ref.vq_flx{2*d};
    W = diag(ref.wq_flx{2*d});
    vals = M \ (P * W * F) / hx(d);
    vals = vals(:);
    res{2*d} = sparse(rows, cols, vals, ng, 1);
end

end


