function res = assemble_advection(par, ref, x, hx)
%ASSEMBLE_ADVECTION
%   Detailed explanation goes here

%========================================================================
% Parse parameters
%========================================================================
ng = ref.n_dofs * prod(par.nx);

%========================================================================
% Assemble functions
%========================================================================
assemble_fn = cell(1, 2);
assemble_fn{1} = { ...
    @assemble_advection_neg, ...
    @assemble_advection_pos, ...
    @assemble_advection_ctr, ...    
};
assemble_fn{2} = @assemble_advection_sl;

%========================================================================
% Advection matrices
%========================================================================
if par.adv_t == 1
    if ismember(par.adv_flx, [1, 2, 3])
        A = cell(1, par.dim);
        for d = 1:par.dim
            A{d} = assemble_fn{1}{par.adv_flx(d)}(par, ref, d, hx);
        end
        res = sparse(ng, ng);
        for d = 1:par.dim
            res = res + par.advection(d) * A{d};
        end
    end
else
    res = cell(1, par.ord_t);
    for i = 1:par.ord_t
        res{i} = assemble_fn{2}(par, ref, i, 1);
    end
end


end

