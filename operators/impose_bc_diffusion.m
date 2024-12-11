function res = impose_bc_diffusion(par, ref, x, hx, t, vec_bdr, mat_aux)
%IMPOSE_BC_DIFFUSION
%   Impose boundary condition for diffusion terms.

%========================================================================
% Parse parameters and initialize variables
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = ref.n_dofs; % Number of local DOFs per cell
ng = nl * nc; % Total number of global DOFs

%========================================================================
% Combine boundary values
%========================================================================
res = zeros(ng, 1);
for d1 = 1:par.dim
    T = sparse(ng, 1);
    for d2 = 1:par.dim
        T = T + par.diffusion(d1, d2) * (vec_bdr{2*d2} - vec_bdr{2*d2-1});
    end
    res = res + mat_aux{d1} * T;

    for d2 = 1:par.dim
        d = par.diffusion(d1, d2);
        h = hx(d2);
        if par.dfn_flx1 == 1
            res = res + d / h * vec_bdr{2*d2-1};
        elseif par.dfn_flx1 == 2
            res = res + d / h * vec_bdr{2*d2};
        elseif par.dfn_flx1 == 3
            res = res + d / h * (vec_bdr{2*d2-1} + ...
                vec_bdr{2*d2}) / 2;
        end
    end
end

end
