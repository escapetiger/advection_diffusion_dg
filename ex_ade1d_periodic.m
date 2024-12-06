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
    'name', '1D Accuracy Test', ... % name of example
    'ax', [-L, L], ... % spatial domain
    'fn_ic', @fn_ic, ... % initial condition
    'fn_exact', @fn_exact, ... % exact solution
    'advection', 0.3, ... % advection coefficients
    'diffusion', 0, ... % diffusion coefficients
    'wavelen', pi/L, ... % wave length
    'amplitude', [0, 1], ... % amplitude
    'bc', 0, ... % boundary condition
    'nx', {8, 16, 32, 64, 128}, ... % number of grid cells in each dimension
    'cfl', 0.5, ... % CFL number
    'ord_t', 2, ... % temporal order
    'ord_x', 2, ... % spatial order
    'poly_t', 'P', ... % polynomial type
    'basis_t', 1, ... % basis type
    'adv_flx', 1, ... % advection flux type
    'dfn_flx1', 3, ... % diffusion flux type for auxiliary variable
    'dfn_flx2', 3, ... % diffusion flux type for primal variable
    't_plot', linspace(0, L, 40), ... % timepoints to plot
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
    fprintf("\n[Simulate with nx = %d]\n", prob(i).nx);
    par{i} = ade_setup(prob(i));
    res{i} = ade_solver(par{i});
    errors{1, i} = res{i}.error.L1;
    errors{2, i} = res{i}.error.L2;
    errors{3, i} = res{i}.error.Linf;
    resolutions(i) = norm((prob(i).ax(2:2:end) - prob(i).ax(1:2:end))./prob(i).nx);
    descriptions{i} = sprintf('%d', prob(i).nx);
    order = min([order, prob(i).ord_t, prob(i).ord_x]);
end
error_table = convergence( ...
    errors, resolutions, descriptions, {'L1', 'L2', 'Linf'}, ones(1, 3)*order);
disp(error_table);

%========================================================================
% Problem Specific Functions
%========================================================================
function f = fn_ic(par, x)
f = par.amplitude(1) * cos(par.wavelen*x{1}) + ...
    par.amplitude(2) * sin(par.wavelen*x{1});
end

function f = fn_exact(par, x, t)
f = exp(-par.lambda*t) * ( ...
    par.amplitude(1) * cos(par.wavelen*(x{1}-par.advection*t)) + ...
    par.amplitude(2) * sin(par.wavelen*(x{1}-par.advection*t)));
end

function output(par, x, U, step)
clf;
Ue = par.fn_exact(par, x, par.t_plot(step));
hold on; plot(x{1},U,'bo'); plot(x{1},Ue,'r-');
xlim(par.ax(1:2)), xlabel('x'), ylim([-1, 1]);
ylabel('u');
hold off;
title(sprintf('%s with DG-%d-%d at t = %0.2f', ...
    par.name, par.ord_t, par.ord_x, par.t_plot(step)));
drawnow;
end
