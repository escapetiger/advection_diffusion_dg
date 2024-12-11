function res = assemble_bc_sl(par, ref, x, hx, t, ht)
%ASSEMBLE_BC_SL
%   Assemble vectors for boundary conditions for semi-Lagrangian terms.
%   Computational domain must be rectangular.

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
nq = ref.nq_vol; % Number of flux quadrature points
M = ref.v_u_vol; % Mass Matrix

res = cell(1, ref.ns);
for i = 1:ref.ns
    res{i} = zeros(ng, 1);
end

for i = 1:ref.ns
    for p = 1:ref.n_volps(i)
        % Find all the elements whose pth volume piece goes outside of
        % computational region.
        m = multi_index(par.nx);
        s = ref.volp_shift{i}(:, p)';
        bi = (m + s < 1) .* (2 * (1:par.dim) - 1) + (m + s > par.nx) .* (2 * (1:par.dim));
        b = m(sum(bi~=0, 2) > 0, :);
        nb = size(b, 1);
        k = m2i(b, par.nx);
        rows = repmat((k - 1)'*nl, nl, 1) + repmat((1:nl)', 1, nb);
        rows = rows(:);
        cols = ones(size(rows));

        if nb == 0, continue; end

        % Find all the quadrature points to evaluate.
        xq = squeeze(ref.xq_volp_e{i}(:, :, p));
        X1 = cell(1, par.dim);
        X2 = cell(1, par.dim);
        for d = 1:par.dim
            X2{d} = reshape(x{d}(b(:, d)), 1, []) + xq(d, :)' * hx(d);
            X2{d} = X2{d}(:);
            X1{d} = X2{d}(:) - ref.dx{i}(d);
        end

        % Find intersection of trajectories and boundaries.
        Xi = cell(1, par.dim); % space coordinates of valid intersection
        ti = cell(1, par.dim); % time coordinates of valid intersection
        fi = cell(1, par.dim); % flags indicates invalid intersection
        ai = par.ax(2*(1:par.dim)-(par.advection>0));
        for d = 1:par.dim
            [Xi{d}, ti{d}, fi{d}] = line_slice_intersect(X1, 0, X2, ref.dt{i}, ai(d), d);
        end
        Xb = cell(1, par.dim);
        tb = zeros(1, nb*nq);
        kb = zeros(1, nb*nq);
        for d = 1:par.dim
            Xb{d} = zeros(1, nb*nq);
            for dd = 1:par.dim
                Xb{d}(~fi{dd}) = Xi{dd}{d};
                tb(~fi{dd}) = ht - (ref.dt{i} - ti{dd});
                kb(~fi{dd}) = ai(dd);
            end
        end

        % Evaluation
        U = par.fn_bc(par, Xb, t + tb, kb);
        U = reshape(U, nq, nb);
        P = squeeze(ref.vq_volp_e{i}(:,:,p));
        W = diag(ref.wq_volp_e{i}(:,p));
        vals = M \ (P * W * U);
        vals = vals(:);
        res{i} = res{i} + sparse(rows, cols, vals, ng, 1);
    end
end


end
