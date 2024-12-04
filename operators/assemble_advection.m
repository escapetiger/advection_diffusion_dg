function res = assemble_advection(par, dg, x, hx)
%ASSEMBLE_ADVECTION
%   Assemble sparse matrices for advection terms.

%========================================================================
% Parse parameters
%========================================================================
ng = dg.n_dofs * prod(par.nx);

%========================================================================
% Assemble functions
%========================================================================
assemble_fn = { ...
    @assemble_advection_neg, ...
    @assemble_advection_pos, ...
    @assemble_advection_ctr};

%========================================================================
% Advection matrices
%========================================================================
A = cell(1, par.dim);
for d = 1:par.dim
    A{d} = assemble_fn{par.adv_flx(d)}(par, dg, d, hx);
end

res = sparse(ng, ng);
for d = 1:par.dim
    res = res + par.advection(d) * A{d};
end


end

