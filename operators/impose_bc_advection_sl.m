function res = impose_bc_advection_sl(par, ref, x, hx, t, vec_bdr, kt)
%IMPOSE_BC_ADVECTION_SL
%   Impose boundary condition for SL advection terms.

%========================================================================
% Parse parameters and initialize variables
%========================================================================
nc = prod(par.nx); % Total number of cells
nl = ref.n_dofs; % Number of local DOFs per cell
ng = nl * nc; % Total number of global DOFs
nv = ref.nv;
nt = ref.nt;

%========================================================================
% Combine boundary values
%========================================================================
res = zeros(ng, 1);
for i = 1:nv
    res = res + vec_bdr{m2i([kt,i], [nt,nv])};
end


end
