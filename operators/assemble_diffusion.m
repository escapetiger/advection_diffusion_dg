function D = assemble_diffusion(par, dg, x, hx)
%ASSEMBLE_DIFFUSION
%   Assemble sparse matrix for diffusion terms using local DG method.

%========================================================================
% Parse parameters
%========================================================================
ng = dg.n_dofs * prod(par.nx);

%========================================================================
% Assemble functions
%========================================================================
assemble_aux_fn = { ...
    @assemble_diffusion_aux_neg, ...
    @assemble_diffusion_aux_pos, ...
    @assemble_diffusion_aux_ctr};

assemble_pri_fn = { ...
    @assemble_diffusion_pri_neg, ...
    @assemble_diffusion_pri_pos, ...
    @assemble_diffusion_pri_ctr};

%========================================================================
% Diffusion matrices
%========================================================================
A1 = cell(1, par.dim);
for d = 1:par.dim
    A1{d} = assemble_aux_fn{par.dfn_flx1}(par, dg, d, hx);
end


A2 = cell(1, par.dim);
for d = 1:par.dim
    A2{d} = assemble_pri_fn{par.dfn_flx2}(par, dg, d, hx);
end

D = sparse(ng, ng);
for d1 = 1:par.dim
    for d2 = 1:par.dim
        if par.diffusion(d1, d2) ~= 0
            D = D + par.diffusion(d1, d2) * A1{d1} * A2{d2};
        end
    end
end

end
