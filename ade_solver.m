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
    nt = ceil(par.tfinal/(par.cfl * max(hx.^(3 / 2))));
else
    nt = ceil(par.tfinal/(par.cfl * max(hx)));
end
if nt > 0
    ht = par.tfinal / nt;
else
    ht = 0;
end

%========================================================================
% Static data (space only)
%========================================================================
dg = setup_dg(par.dim, par.ord_x, par.basis_t, par.poly_t, par.n_plot, ...
    par.n_error);

% identity matrix
ng = dg.n_dofs * prod(par.nx);
mat_eye = speye(ng, ng);
% advection matrix
tic
fprintf("Assembling advection matrix...\n");
mat_adv = assemble_advection(par, dg, x, hx);
cputime(1) = cputime(1) + toc;

% diffusion matrix
tic
fprintf("Assembling diffusion matrix...\n");
mat_dfn = assemble_diffusion(par, dg, x, hx);

switch par.ord_t
    case 1
        A = (mat_eye - ht * (mat_dfn - mat_adv));
    case 2
        A1 = (mat_eye - ht * (mat_dfn - mat_adv));
        A2 = (3 * mat_eye - 2 * ht * (mat_dfn - mat_adv));
    case 3
        A1 = (mat_eye - ht * (mat_dfn - mat_adv));
        A2 = (3 * mat_eye - 2 * ht * (mat_dfn - mat_adv));
        A3 = (11 * mat_eye - 6 * ht * (mat_dfn - mat_adv));
end
is_spd = false;
if par.ord_x == 1
    if issymmetric(A)
        [~, p] = chol(A);
        is_spd = (p == 0);
    end
end
cputime(1) = cputime(1) + toc;

%========================================================================
% Time loop
%========================================================================
par.t_plot = [par.t_plot, inf];
t = 0;
step_count = 0;
plot_count = 1;

% impose initial condition
tic
U = impose_ic(par, dg, x, hx);
cputime(1) = cputime(1) + toc;

fprintf("Simulating...\n");
while t < par.tfinal && step_count < nt
    % update g
    if any(par.bc > 0)
        g = assemble_bc(par, dg, x, hx, t, ht);
    end

    % update U
    tic;
    switch par.ord_t
        case 1
            U0 = U;
            U = spsolve(A, U0, struct('is_spd', is_spd));
        case 2
            if step_count == 0
                U0 = U;
                U = spsolve(A1, U0, struct('is_spd', is_spd));
                U1 = U;
            else                
                U = spsolve(A2, -U0 + 4 * U1, struct('is_spd', is_spd));
                U0 = U1;
                U1 = U;
            end
        case 3
            % NOTE: The first two steps are predicted by BE and BDF2.
            % However, this treatment leads to order reduction in BDF3
            % when using hyperbolic time step.
            if step_count == 0
                U0 = U;
                U = spsolve(A1, U0, struct('is_spd', is_spd));
                U1 = U;
            elseif step_count == 1
                U = spsolve(A2, -U0 + 4 * U1, struct('is_spd', is_spd));
                U2 = U;
            else
                U = spsolve(A3, 2 * U0 - 9 * U1 + 18 * U2, struct('is_spd', is_spd));
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
                    U_plot0 = evaluate(par, dg, U0(:, par.outvar));
                    U_plot = evaluate(par, dg, U(:, par.outvar));
                    U_plot = -lambda * U_plot0 + (1 + lambda) * U_plot;
                case 2
                    if step_count == 0
                        U_plot0 = evaluate(par, dg, U0(:, par.outvar));
                        U_plot = evaluate(par, dg, U(:, par.outvar));
                        U_plot = -lambda * U_plot0 + (1 + lambda) * U_plot;
                    else
                        U_plot0 = evaluate(par, dg, U0(:, par.outvar));
                        U_plot1 = evaluate(par, dg, U1(:, par.outvar));
                        U_plot = evaluate(par, dg, U(:, par.outvar));
                        U_plot = (1 + lambda) * lambda / 2 * U_plot0 - ...
                            lambda * (2 + lambda) * U_plot1 + ...
                            (1 + lambda) * (1 + lambda / 2) * U_plot;
                    end
                case 3
                    if step_count == 0
                        U_plot0 = evaluate(par, dg, U0(:, par.outvar));
                        U_plot = evaluate(par, dg, U(:, par.outvar));
                        U_plot = -lambda * U_plot0 + (1 + lambda) * U_plot;
                    elseif step_count == 1
                        U_plot0 = evaluate(par, dg, U0(:, par.outvar));
                        U_plot1 = evaluate(par, dg, U1(:, par.outvar));
                        U_plot = evaluate(par, dg, U(:, par.outvar));
                        U_plot = (1 + lambda) * lambda / 2 * U_plot0 - ...
                            lambda * (2 + lambda) * U_plot1 + ...
                            (1 + lambda) * (1 + lambda / 2) * U_plot;
                    else
                        U_plot0 = evaluate(par, dg, U0(:, par.outvar));
                        U_plot1 = evaluate(par, dg, U1(:, par.outvar));
                        U_plot2 = evaluate(par, dg, U2(:, par.outvar));
                        U_plot = evaluate(par, dg, U(:, par.outvar));
                        U_plot = -(lambda + 2) * (lambda + 1) * lambda / 6 * U_plot0 + ...
                            (lambda + 3) * (lambda + 1) * lambda / 2 * U_plot1 - ...
                            (lambda + 3) * (lambda + 2) * lambda / 2 * U_plot2 + ...
                            (lambda + 3) * (lambda + 2) * (lambda + 1) / 6 * U_plot;
                    end
            end
            x_plot = cell(1, par.dim);
            for d = 1:par.dim
                x_plot{d} = repmat(x{d}, dg.np, 1) + repmat(dg.xp{d}', 1, par.nx(d)) * hx(d);
                x_plot{d} = x_plot{d}(:);
            end
            if par.dim > 1
                U_plot = reshape(U_plot, par.nx.*dg.np);
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
    error = evaluate_error(par, dg, U(:, par.outvar), x, hx, t);
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
