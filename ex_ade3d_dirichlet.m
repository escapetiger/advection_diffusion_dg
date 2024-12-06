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
    'name', '3D Accuracy Test', ... % name of example
    'ax', [-L, L, -L, L, -L, L], ... % spatial domain
    'fn_ic', @fn_ic, ... % initial condition
    'fn_bc', @fn_bc, ... % boundary condition
    'fn_exact', @fn_exact, ... % exact solution
    'advection', [0.5, 0.5, 0.5], ... % advection coefficients
    'diffusion', [0.05, -0.01, 0.01; 0.01, 0.06, 0.01;0.01, 0.01, 0.07], ... % diffusion coefficients
    'wavelen', [pi / L, pi / L, pi / L], ... % wave length
    'bc', [1, 1, 1], ... % boundary condition
    'nx', {[4, 3, 2], [8, 6, 4], [16, 12, 8], [32, 24, 16]}, ... % number of grid cells in each dimension
    'cfl', 0.5, ... % CFL number
    'ord_t', 3, ... % temporal order
    'ord_x', 3, ... % spatial order
    'poly_t', 'P', ... % polynomial type
    'basis_t', 1, ... % basis type
    'adv_flx', [1, 1, 1], ... % advection flux type
    'dfn_flx1', 2, ... % diffusion flux type for auxiliary variable
    'dfn_flx2', 1, ... % diffusion flux type for primal variable
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
    fprintf("\n[Simulate with nx = (%d, %d, %d)]\n", prob(i).nx);
    par{i} = ade_setup(prob(i));
    res{i} = ade_solver(par{i});
    errors{1, i} = res{i}.error.L1;
    errors{2, i} = res{i}.error.L2;
    errors{3, i} = res{i}.error.Linf;
    resolutions(i) = norm((prob(i).ax(2:2:end) - prob(i).ax(1:2:end))./prob(i).nx);
    descriptions{i} = sprintf('%dx%dx%d', prob(i).nx);
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
y1 = x{1} - par.advection(1) * t;
y2 = x{2} - par.advection(2) * t;
y3 = x{3} - par.advection(3) * t;
z1 = par.wavelen(1) * y1;
z2 = par.wavelen(2) * y2;
z3 = par.wavelen(3) * y3;
f = exp(-par.lambda*t) * sin(z1 + z2 + z3);
end


function output(par, x, U, step)
% Plotting routine, showing the 3D results and a 1D cross-section.
t = par.t_plot(step); % Time of output]7
X = cell(1, par.dim);
[X{:}] = ndgrid(x{:});
Ue = par.fn_exact(par, X, t);
U = permute(U, [2, 1, 3]); 
Ue = permute(Ue, [2, 1, 3]);
xslice = 0;
yslice = 0;
zslice = 0; % Slice locations.
clf, subplot(1, 3, 1:2);
h = slice(x{1}, x{2}, x{3}, U, xslice, yslice, zslice);
set(h, 'EdgeColor', 'none');
xlabel('x'), ylabel('y'), zlabel('z');
axis vis3d equal, caxis([-1 1]);
title(sprintf('%s with DG-%d-%d at t = %0.2f', par.name, par.ord_t, par.ord_x, t));
colormap jet(255); colorbar;

subplot(1, 3, 3);
r = linspace(min(x{1}), max(x{1}), 64);
% 1D cross sections
Uexyz = interp3(x{1}, x{2}, x{3}, Ue, r/sqrt(3), r/sqrt(3), r/sqrt(3));
Uxyz = interp3(x{1}, x{2}, x{3}, U, r/sqrt(3), r/sqrt(3), r/sqrt(3));
% Plot 1D cross sections.
plot(r, Uexyz, 're', 'DisplayName', 'Ref');
hold on;
plot(r, Uxyz, 'bo', 'DisplayName', 'DG');
hold off;
axis([min(x{1}), max(x{1}), [-1, 1]]);
legend('Location', 'SouthEast');
title('Cut at x = y = z');
drawnow
end
