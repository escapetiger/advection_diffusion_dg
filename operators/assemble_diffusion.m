function [D, A1, A2] = assemble_diffusion(par, ref, x, hx)
%ASSEMBLE_DIFFUSION
%   Assemble linear diffusion matrix using local DG method.

%========================================================================
% Parse parameters
%========================================================================
ng = ref.n_dofs * prod(par.nx);

%========================================================================
% Assemble functions
%========================================================================
assemble_aux_fn = { ...
    @assemble_diffusion_aux_neg, ...
    @assemble_diffusion_aux_pos, ...
    @assemble_diffusion_aux_ctr};

assemble_aux_bdr_ext_fn = { ...
    @assemble_diffusion_aux_bdr_ext_neg, ...
    @assemble_diffusion_aux_bdr_ext_pos, ...
    @assemble_diffusion_aux_bdr_ext_ctr};

assemble_aux_bdr_jmp_fn = { ...
    @assemble_diffusion_aux_bdr_jmp_neg, ...
    @assemble_diffusion_aux_bdr_jmp_pos, ...
    @assemble_diffusion_aux_bdr_jmp_ctr};

assemble_pri_fn = { ...
    @assemble_diffusion_pri_neg, ...
    @assemble_diffusion_pri_pos, ...
    @assemble_diffusion_pri_ctr};

%========================================================================
% Diffusion matrices
%========================================================================
A1 = cell(1, par.dim);
for d = 1:par.dim
    A1{d} = assemble_aux_fn{par.dfn_flx1}(par, ref, d, hx);
    if par.bc(d) == 1
        A1_ext = assemble_aux_bdr_ext_fn{par.dfn_flx1}(par, ref, d, hx);
        A1{d} = A1{d} + A1_ext;
    end
end

B1 = cell(1, par.dim);
for d = 1:par.dim
    if par.bc(d) == 1
        B1{d} = assemble_aux_bdr_jmp_fn{par.dfn_flx1}(par, ref, d, hx);
    end
end

A2 = cell(1, par.dim);
for d = 1:par.dim
    A2{d} = assemble_pri_fn{par.dfn_flx2}(par, ref, d, hx);
end

D = sparse(ng, ng);
for d1 = 1:par.dim
    for d2 = 1:par.dim
        if par.diffusion(d1, d2) ~= 0
            D = D + par.diffusion(d1, d2) * A1{d1} * A2{d2};
        end
    end
end

for d1 = 1:par.dim
    if par.bc(d1) == 1
        for d2 = 1:par.dim
            d = par.diffusion(d1, d2);
            h = hx(d2);
            if par.dfn_flx1 == 1 && par.dfn_flx2 == 2
                D = D + d / h * B1{d2};
            elseif par.dfn_flx1 == 2 && par.dfn_flx2 == 1
                D = D + d / h * B1{d2};
            elseif par.dfn_flx1 == 3 && par.dfn_flx2 == 3
                D = D + d / h * B1{d2};
            end
        end
    end
end

if all(par.bc == 0)
    A1 = [];
    A2 = [];
end

end
