function res = impose_bc_advection(par, dg, x, hx, t, vec_bdr)
%IMPOSE_BC_ADVECTION
%   Impose boundary condition for advection terms.

%========================================================================
% Parse parameters and initialize variables
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = dg.n_dofs; % Number of local DOFs per cell
ng = nl * nc; % Total number of global DOFs

%========================================================================
% Combine boundary values
%========================================================================
res = zeros(ng, 1);
for d = 1:par.dim
    switch par.adv_flx(d)
        case 1
            res = res - par.advection(d) * vec_bdr{2*d-1};
        case 2
            res = res + par.advection(d) * vec_bdr{2*d};
    end
end

end
