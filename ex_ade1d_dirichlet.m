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
    'name', '1D Dirichlet Test', ... % name of example
    'ax', [0, L], ... % spatial domain
    'fn_ic', @fn_ic, ... % initial condition
    'fn_bc', @fn_bc, ... % boundary condition
    'fn_exact', @fn_exact, ... % exact solution
    'advection', 0.5, ... % advection coefficients
    'diffusion', 0.1, ... % diffusion coefficients
    'wavelen', pi/L, ... % wave length
    'amplitude', 1, ... % amplitude
    'bc', 1, ... % boundary condition
    'nx', {8, 16, 32, 64, 128, 256}, ... % number of grid cells in each dimension
    'cfl', 0.3, ... % CFL number
    'ord_t', 2, ... % temporal order
    'ord_x', 2, ... % spatial order
    'poly_t', 'P', ... % polynomial type
    'basis_t', 1, ... % basis type
    'adv_flx', 1, ... % advection flux type
    'dfn_flx1', 2, ... % diffusion flux type for auxiliary variable
    'dfn_flx2', 1, ... % diffusion flux type for prime variable
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
f = fn_exact(par, x, 0);
end

function f = fn_bc(par, x, t)
f = fn_exact(par, x, t);
end

function f = fn_exact(par, x, t)
y = x{1} - par.advection * t;
z = par.wavelen * y;
f = exp(-par.lambda*t) * par.amplitude * sin(z);
end

function output(par, x, U, step)
clf;
hold on; plot(x{1},U,'bo'); 
Ue = par.fn_exact(par, x, par.t_plot(step));
plot(x{1},Ue,'r-');
xlim(par.ax(1:2)), xlabel('x');
ylabel('u');
hold off;
title(sprintf('%s with DG-%d-%d at t = %0.2f', ...
    par.name, par.ord_t, par.ord_x, par.t_plot(step)));
drawnow;
end
