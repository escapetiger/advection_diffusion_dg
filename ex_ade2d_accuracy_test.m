%========================================================================
% Environment
%========================================================================
clc;
clear;
addpath('utilities');
addpath('operators');

%========================================================================
% Problem Parameters
%========================================================================
L = 1;
prob = struct( ...
    'name', '2D Accuracy Test', ... % name of example
    'ax', [-L, L, -L, L], ... % spatial domain
    'fn_ic', @fn_ic, ... % initial condition
    'fn_exact', @fn_exact, ... % exact solution
    'advection', [0.5, 0.5], ... % advection coefficients
    'diffusion', [0.05, -0.02; -0.02, 0.05], ... % diffusion coefficients
    'wavelen', [pi / L, pi / L], ... % wave length
    'amplitude', [0, 1], ... % amplitude
    'bc', [0, 0], ... % boundary condition
    'nx', {[8, 8], [16, 16], [32, 32], [64, 64]}, ... % number of grid cells in each dimension
    'cfl', 0.5, ... % CFL number
    'ord_t', 2, ... % temporal order
    'ord_x', 2, ... % spatial order
    'poly_t', 'Q', ... % polynomial type
    'basis_t', 2, ... % basis type
    'adv_flx', [1, 1], ... % advection flux type
    'dfn_flx1', 1, ... % diffusion flux type for auxiliary variable
    'dfn_flx2', 2, ... % diffusion flux type for primal variable
    't_plot', linspace(0, L, 20), ... % timepoints to plot
    'output', @output, ... % customized output routine
    'outvar', 1, ... % output variable
    'n_plot', 1, ... % number of spacepoints to plot
    'n_error', 10 ... % number of points used in error computation
    );

%========================================================================
% Simulation Execution
%========================================================================
par = cell(1, numel(prob));
res = cell(1, numel(prob));
errors = cell(3, numel(prob));
resolutions = ones(1, numel(prob));
descriptions = cell(1, numel(prob));
order = inf;
for i = 1:numel(prob)
    fprintf("\n[Simulate with nx = (%d, %d)]\n", prob(i).nx);
    par{i} = ade_setup(prob(i));
    res{i} = ade_solver(par{i});
    errors{1, i} = res{i}.error.L1;
    errors{2, i} = res{i}.error.L2;
    errors{3, i} = res{i}.error.Linf;
    resolutions(i) = norm((prob(i).ax(2:2:end) - prob(i).ax(1:2:end))./prob(i).nx);
    descriptions{i} = sprintf('%dx%d', prob(i).nx);
    order = min([order, prob(i).ord_t, prob(i).ord_x]);
end
error_table = convergence( ...
    errors, resolutions, descriptions, {'L1', 'L2', 'Linf'}, ones(1, 3)*order);
disp(error_table);

%========================================================================
% Problem Specific Functions
%========================================================================
function f = fn_ic(par, x)
z = par.wavelen(1) * x{1} + par.wavelen(2) * x{2};
f = par.amplitude(1) * cos(z) + par.amplitude(2) * sin(z);
end

function f = fn_exact(par, x, t)
z = par.wavelen(1) * (x{1} - par.advection(1) * t) + ...
    par.wavelen(2) * (x{2} - par.advection(2) * t);
f = exp(-par.lambda*t) * (par.amplitude(1) * cos(z) + ...
    par.amplitude(2) * sin(z));
end

function output(par, x, U, step)
t = par.t_plot(step);

clf, subplot(1, 3, 1:2);
imagesc(x{1}, x{2}, U');
axis xy equal tight;
xlim(par.ax(1:2)), xlabel('x');
ylim(par.ax(3:4)), ylabel('y');
colormap jet(255); colorbar; caxis([-1 1]);
title(sprintf('%s with DG-%d-%d at t = %0.2f', ...
    par.name, par.ord_t, par.ord_x, t));

subplot(1, 3, 3);
X = cell(1, par.dim);
[X{:}] = ndgrid(x{:});
Ue = par.fn_exact(par, X, t);
r = linspace(min(x{1}), max(x{1}), 64);
Uexy = interp2(x{1}, x{2}, Ue, r/sqrt(2), r/sqrt(2));
Uxy = interp2(x{1}, x{2}, U, r/sqrt(2), r/sqrt(2));
plot(r, Uexy, 're', 'DisplayName', 'Ref');
hold on;
plot(r, Uxy, 'bo', 'DisplayName', 'DG');
hold off;
axis([min(x{1}), max(x{1}), [-1, 1]]);
legend('Location', 'SouthEast');
title('Cut at x = y');
drawnow;
end