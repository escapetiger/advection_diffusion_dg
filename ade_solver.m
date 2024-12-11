function res = ade_solver(par)
%DFN_SOLVER
%   Linear diffusion solver.

%========================================================================
% Global monitor
%========================================================================
cputime = zeros(1, 3);

%========================================================================
% Geometric discretization
%========================================================================

% spatial grid
hx = (par.ax(2:2:end) - par.ax(1:2:end)) ./ par.nx;
switch par.dim
    case 1
        x{1} = par.ax(1) + hx(1) / 2:hx(1):par.ax(2) - hx(1) / 2;
    case 2
        x{1} = par.ax(1) + hx(1) / 2:hx(1):par.ax(2) - hx(1) / 2;
        x{2} = par.ax(3) + hx(2) / 2:hx(2):par.ax(4) - hx(2) / 2;
    case 3
        x{1} = par.ax(1) + hx(1) / 2:hx(1):par.ax(2) - hx(1) / 2;
        x{2} = par.ax(3) + hx(2) / 2:hx(2):par.ax(4) - hx(2) / 2;
        x{3} = par.ax(5) + hx(3) / 2:hx(3):par.ax(6) - hx(3) / 2;
end

% temporal grid
if par.ord_t == 3
    % NOTE: In BDF3, the first two steps are predicted by BE and BDF2.
    % However, this treatment leads to order reduction when using 
    % hyperbolic time step ht=C*hx. We have to use time step 
    % ht=C*hx^(3/2) to get desired convergence order.
    nt = ceil(par.tfinal/(par.cfl * min(hx.^(3/2))));
else
    nt = ceil(par.tfinal/(par.cfl * min(hx)));
end
if nt > 0
    ht = par.tfinal / nt;
else
    ht = 0;
end

%========================================================================
% Static data
%========================================================================
ref = struct;
ref = setup_dg(ref, par.dim, par.ord_x, par.basis_t, par.poly_t, ...
    par.n_plot, par.n_error);
if par.adv_t == 2
    ref = setup_sldg(ref, par.dim, par.ord_x, par.basis_t, ...
        par.poly_t, par.advection, ht*(1:par.ord_t), hx);
end

% Advection operator
tic
fprintf("Assembling advection operator...\n");
mat_adv = assemble_advection(par, ref, x, hx);
cputime(1) = cputime(1) + toc;


% Diffusion operator
tic
fprintf("Assembling diffusion operator...\n");
[mat_dfn, mat_dfn_aux, ~] = assemble_diffusion(par, ref, x, hx);
cputime(1) = cputime(1) + toc;

% Evolution operator
tic
fprintf("Assembling evolution operator...\n");
[A, is_A_spd] = assemble_evolution(par, ref, ht, mat_adv, mat_dfn);
if par.adv_t == 1
    mat_adv = [];
end
cputime(1) = cputime(1) + toc;


%========================================================================
% Time loop
%========================================================================
par.t_plot = [par.t_plot, inf];
t = 0;
step_count = 0;
plot_count = 1;

% Impose initial condition
tic
U = impose_ic(par, ref, x, hx);
cputime(1) = cputime(1) + toc;

fprintf("Simulating...\n");
while t < par.tfinal && step_count < nt
    % Update boundary condition: Ub
    Ub = assemble_bc(par, ref, x, hx, t, ht);
    if par.adv_t == 2
        Ub_sl = assemble_bc_sl(par, ref, x, hx, t, ht);
    else
        Ub_sl = [];
    end

    % Update solution: U
    tic;
    switch par.ord_t
        case 1
            U0 = U;
            U = step_BE( ...
                par, ref, x, hx, t, ht, U0, A{1}, is_A_spd, ...
                Ub, Ub_sl, mat_adv, mat_dfn_aux);
        case 2
            if step_count == 0
                U0 = U;
                U = step_BE( ...
                    par, ref, x, hx, t, ht, U0, A{1}, is_A_spd{1}, ...
                    Ub, Ub_sl, mat_adv, mat_dfn_aux);
                U1 = U;
            else                
                U = step_BDF2( ...
                    par, ref, x, hx, t, ht, U0, U1, A{2}, is_A_spd{2}, ...
                    Ub, Ub_sl, mat_adv, mat_dfn_aux);
                U0 = U1;
                U1 = U;
            end
        case 3
            if step_count == 0
                U0 = U;
                U = step_BE( ...
                    par, ref, x, hx, t, ht, U0, A{1}, is_A_spd{1}, ...
                    Ub, Ub_sl, mat_adv, mat_dfn_aux);
                U1 = U;
            elseif step_count == 1
                U = step_BDF2( ...
                    par, ref, x, hx, t, ht, U0, U1, A{2}, is_A_spd{2}, ...
                    Ub, Ub_sl, mat_adv, mat_dfn_aux);
                U2 = U;
            else
                U = step_BDF3( ...
                    par, ref, x, hx, t, ht, U0, U1, U2, A{3}, is_A_spd{3}, ...
                    Ub, Ub_sl, mat_adv, mat_dfn_aux);
                U0 = U1;
                U1 = U2;
                U2 = U;
            end
    end
    cputime(1) = cputime(1) + toc;

    % time increment
    t = t + ht;

    % plotting
    tic;
    if par.n_plot > 0
        while t >= par.t_plot(plot_count) - 1e-14
            lambda = (par.t_plot(plot_count) - t) / ht;
            switch par.ord_t
                case 1
                    U_plot = plot_BE(par, ref, lambda, U0, U);
                case 2
                    if step_count == 0
                        U_plot = plot_BE(par, ref, lambda, U0, U);
                    else
                        U_plot = plot_BDF2(par, ref, lambda, U0, U1, U);
                    end
                case 3
                    if step_count == 0
                        U_plot = plot_BE(par, ref, lambda, U0, U);
                    elseif step_count == 1
                        U_plot = plot_BDF2(par, ref, lambda, U0, U1, U);
                    else
                        U_plot = plot_BDF3(par, ref, lambda, U0, U1, U2, U);
                    end
            end
            x_plot = cell(1, par.dim);
            for d = 1:par.dim
                x_plot{d} = repmat(x{d}, ref.np, 1) + repmat(ref.xp{d}', 1, par.nx(d)) * hx(d);
                x_plot{d} = x_plot{d}(:);
            end

            if par.dim > 1
                U_plot = reshape(U_plot, par.nx.*ref.np);
            end
            if nargout(par.output)
                par = par.output(par, x_plot, U_plot, plot_count);
            else
                par.output(par, x_plot, U_plot, plot_count)
            end
            plot_count = plot_count + 1;
        end
    end
    cputime(2) = cputime(2) + toc;

    % Finalize time step
    step_count = step_count + 1;

    fprintf("step = %d; time = %f; step size = %f; mesh size = %f;\n", ...
        step_count, t, ht, max(hx));
end

%========================================================================
% Residual computation
%========================================================================
res = struct;
res.x = x;
res.t = t;
res.U = U;

if par.n_error > 0
    tic
    error = evaluate_error(par, ref, U(:, par.outvar), x, hx, t);
    cputime(3) = cputime(3) + toc;
    fprintf(['Residuals:\n', ...
        ' L1:   %8.3e\n', ...
        ' L2:   %8.3e\n', ...
        ' Linf:   %8.3e\n'], ...
        error.L1, error.L2, error.Linf);
    res.error = error;
end

cputime = reshape([cputime; cputime / sum(cputime) * 1e2], 1, []);
fprintf(['CPU-times\n', ...
    ' U-step:%12.2fs%5.0f%%\n', ...
    ' Plotting:%12.2fs%5.0f%%\n', ...
    ' Residual:%12.2fs%5.0f%%\n'], cputime);
end

function [A, spd] = assemble_evolution(par, ref, ht, mat_adv, mat_dfn)
% identity operator
ng = ref.n_dofs * prod(par.nx);
mat_eye = speye(ng, ng);

% Evolution operator
A = cell(1, par.ord_t);
switch par.ord_t
    case 1
        A{1} = mat_eye - ht * mat_dfn;
        if par.adv_t == 1
            A{1} = A{1} + ht * mat_adv;
        end
    case 2
        A{1} = mat_eye - ht * mat_dfn;
        A{2} = 3 * mat_eye - 2 * ht * mat_dfn;
        if par.adv_t == 1
            A{1} = A{1} + ht * mat_adv;
            A{2} = A{2} + 2 * ht * mat_adv;
        end
    case 3
        A{1} = mat_eye - ht * mat_dfn;
        A{2} = 3 * mat_eye - 2 * ht * mat_dfn;
        A{3} = 11 * mat_eye - 6 * ht * mat_dfn;
        if par.adv_t == 1
            A{1} = A{1} + ht * mat_adv;
            A{2} = A{2} + 2 * ht * mat_adv;
            A{3} = A{3} + 6 * ht * mat_adv;
        end
end

spd = cell(1, par.ord_t);
for i = 1:par.ord_t
    spd{i} = false;
    if par.ord_x == 1
        if issymmetric(A{i})
            [~, p] = chol(A{i});
            spd{i} = (p == 0);
        end
    end
end

end

function U = step_BE( ...
    par, ref, x, hx, t, ht, U0, A, is_A_spd, Ub, Ub_sl, mat_adv, mat_dfn_aux)
% Evolve by backward Euler stepper
    switch par.adv_t
        case 1
            b = U0;
            if any(par.bc > 0) && ~isempty(mat_dfn_aux) && ~isempty(Ub)
                b0 = impose_bc_diffusion(par, ref, x, hx, t, Ub, mat_dfn_aux);
                b1 = impose_bc_advection(par, ref, x, hx, t, Ub);
                b = b + ht * (b0 - b1);
            end
        case 2
            b = mat_adv{1} * U0;
            if any(par.bc > 0) && ~isempty(mat_dfn_aux) && ~isempty(Ub) && ~isempty(Ub_sl)
                b0 = impose_bc_diffusion(par, ref, x, hx, t, Ub, mat_dfn_aux);
                b1 = impose_bc_advection_sl(par, ref, x, hx, t, Ub_sl, 1);
                b = b + ht * b0 + b1;
            end
    end
    U = spsolve(A, b, struct('is_spd', is_A_spd));
end

function U = step_BDF2( ...
    par, ref, x, hx, t, ht, U0, U1, A, is_A_spd, Ub, Ub_sl, mat_adv, mat_dfn_aux)
% Evolve by the second order backward differential formula 
    switch par.adv_t
        case 1
            b = -U0 + 4 * U1;
            if any(par.bc > 0)
                b0 = impose_bc_diffusion(par, ref, x, hx, t, Ub, mat_dfn_aux);
                b1 = impose_bc_advection(par, ref, x, hx, t, Ub);
                b = b + 2 * ht * (b0 - b1);
            end
        case 2
            b = -mat_adv{2} * U0 + 4 * mat_adv{1} * U1;
            if any(par.bc > 0) && ~isempty(mat_dfn_aux) && ~isempty(Ub) && ~isempty(Ub_sl)
                b0 = impose_bc_diffusion(par, ref, x, hx, t, Ub, mat_dfn_aux);
                b1 = impose_bc_advection_sl(par, ref, x, hx, t, Ub_sl, 1);
                b2 = impose_bc_advection_sl(par, ref, x, hx, t, Ub_sl, 2);
                b = b + 2 * ht * b0 + 4 * b1 - b2;
            end
    end
    U = spsolve(A, b, struct('is_spd', is_A_spd));
end

function U = step_BDF3( ...
    par, ref, x, hx, t, ht, U0, U1, U2, A, is_A_spd, Ub, Ub_sl, mat_adv, mat_dfn_aux)
% Evolve by the third order backward differential formula
    switch par.adv_t
        case 1
            b = 2 * U0 - 9 * U1 + 18 * U2;
            if any(par.bc > 0)
                b0 = impose_bc_diffusion(par, ref, x, hx, t, Ub, mat_dfn_aux);
                b1 = impose_bc_advection(par, ref, x, hx, t, Ub);
                b = b + 6 * ht * (b0 - b1);
            end
        case 2
            b = 2 * mat_adv{3} * U0 - 9 * mat_adv{2} * U1 + 18 * mat_adv{1} * U2;
            if any(par.bc > 0) && ~isempty(mat_dfn_aux) && ~isempty(Ub) && ~isempty(Ub_sl)
                b0 = impose_bc_diffusion(par, ref, x, hx, t, Ub, mat_dfn_aux);
                b1 = impose_bc_advection_sl(par, ref, x, hx, t, Ub_sl, 1);
                b2 = impose_bc_advection_sl(par, ref, x, hx, t, Ub_sl, 2);
                b3 = impose_bc_advection_sl(par, ref, x, hx, t, Ub_sl, 3);
                b = b + 6 * ht * b0 + 2 * b3 - 9 * b2 + 18 * b1;
            end
    end
    U = spsolve(A, b, struct('is_spd', is_A_spd));
end

function U_plot = plot_BE(par, ref, lambda, U0, U)
    U_plot0 = evaluate(par, ref, U0(:, par.outvar));
    U_plot = evaluate(par, ref, U(:, par.outvar));
    U_plot = -lambda * U_plot0 + (1 + lambda) * U_plot;
end

function U_plot = plot_BDF2(par, ref, lambda, U0, U1, U)
    U_plot0 = evaluate(par, ref, U0(:, par.outvar));
    U_plot1 = evaluate(par, ref, U1(:, par.outvar));
    U_plot = evaluate(par, ref, U(:, par.outvar));
    U_plot = (1 + lambda) * lambda / 2 * U_plot0 - ...
        lambda * (2 + lambda) * U_plot1 + ...
        (1 + lambda) * (1 + lambda / 2) * U_plot;
end

function U_plot = plot_BDF3(par, ref, lambda, U0, U1, U2, U)
    U_plot0 = evaluate(par, ref, U0(:, par.outvar));
    U_plot1 = evaluate(par, ref, U1(:, par.outvar));
    U_plot2 = evaluate(par, ref, U2(:, par.outvar));
    U_plot = evaluate(par, ref, U(:, par.outvar));
    U_plot = -(lambda + 2) * (lambda + 1) * lambda / 6 * U_plot0 + ...
        (lambda + 3) * (lambda + 1) * lambda / 2 * U_plot1 - ...
        (lambda + 3) * (lambda + 2) * lambda / 2 * U_plot2 + ...
        (lambda + 3) * (lambda + 2) * (lambda + 1) / 6 * U_plot;
end
