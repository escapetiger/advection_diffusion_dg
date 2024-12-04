function par = ade_setup(prob)
%ADE_SETUP
%   Convert user-defined variables to format compatiable with
%   advection-diffusion solvers.

%========================================================================
% Set struct
%========================================================================
par = prob;
par.prob = prob;

%========================================================================
% Set dimension of problem
%========================================================================
% Detect dimension of problem
if isfield(prob, 'ax'), par.dim = length(prob.ax) / 2;
elseif isfield(prob, 'nx'), par.dim = length(prob.nx);
else, error('At least one of ax, nx, or bc should be defined');
end
% Check if variables are compatible
if isfield(prob, 'ax')
    if isfield(prob, 'nx')
        if abs(length(prob.ax)/2-length(prob.nx)) > 0
            error('Incompatible problem variables ax and nx.')
        end
    end
    if isfield(prob, 'bc')
        if abs(length(prob.ax)/2-length(prob.bc)) > 0
            error('Incompatible problem variables ax and bc.')
        end
    end
elseif isfield(prob, 'nx')
    if isfield(prob, 'bc')
        if abs(length(prob.nx)-length(prob.bc)) > 0
            error('Incompatible problem variables nx and bc.')
        end
    end
end

%========================================================================
% Adapt problem-specific data to solvers
%========================================================================
if isfield(prob, 'fn_bc')
    switch nargin(prob.fn_bc)
        case 2
            par.fn_bc = @(par, x, ~, ~) prob.fn_bc(par, x);
        case 3
            par.fn_bc = @(par, x, t, ~) prob.fn_bc(par, x, t);
        case 4
            par.fn_bc = @(par, x, t, i) prob.fn_bc(par, x, t, i);
    end
end

if isfield(prob, 'wavelen')
    if length(prob.wavelen) ~= par.dim
        error('Invalid wave lengths.');
    end
end

if isfield(prob, 'diffusion')
    if any(size(prob.diffusion) ~= par.dim)
        error('Invalid diffusion coefficients.');
    end
end

if isfield(prob, 'wavelen') && isfield(prob, 'diffusion')
    par.wavelen = reshape(prob.wavelen, 1, par.dim);
    par.lambda = par.wavelen * prob.diffusion * par.wavelen';
end

%========================================================================
% Set defaults
%========================================================================
if ~isfield(par, 'name'), par.name = 'Output'; end
if ~isfield(par, 'fn_ic'), par.fn_ic = @zero; end
if ~isfield(par, 'fn_bc'), par.fn_bc = @zero; end
if ~isfield(par, 'bc'), par.bc = [0, 0, 0]; end
if ~isfield(par, 'ax'), par.ax = [0, 1, 0, 1, 0, 1]; end
if ~isfield(par, 'nx'), par.nx = [100, 100, 100]; end
if ~isfield(par, 'cfl'), par.cfl = 0.99; end
if ~isfield(par, 't_plot'), par.t_plot = 0:.1:1; end
if ~isfield(par, 'tfinal'), par.tfinal = par.t_plot(end); end
if ~isfield(par, 'outvar'), par.outvar = 1; end
if ~isfield(par, 'n_plot'), par.n_plot = 1; end
if ~isfield(par, 'n_error'), par.n_error = 0; end
if ~isfield(par, 'fn_exact'), par.fn_exact = @zero; end
if ~isfield(par, 'advection'), par.advection = [1, 1, 1]; end
if ~isfield(par, 'diffusion')
    par.diffusion = [1, 0, 0; 0, 1, 0; 0, 0, 1];
end
switch par.dim
    case 1
        if ~isfield(par, 'output'), par.output = @default_output_1d; end
    case 2
        if ~isfield(par, 'output'), par.output = @default_output_2d; end
    case 3
        if ~isfield(par, 'output'), par.output = @default_output_3d; end
end


end

%========================================================================
% Utility functions
%========================================================================
function f = zero(varargin), f = zeros(size(varargin{1})); end

function default_output_1d(par, x, U, step)
clf, plot(x{1}, U), xlim(par.ax(1:2)), xlabel('x');
title(sprintf('%s with DG-%d-%d at t = %0.2f', ...
    par.name, par.ord_t, par.ord_x, par.t_plot(step)));
drawnow;
end

function default_output_2d(par, x, U, step)
clf, imagesc(x{1}, x{2}, U'), axis xy equal tight;
title(sprintf('%s with DG-%d-%d at t = %0.2f', ...
    par.name, par.ord_t, par.ord_x, par.t_plot(step)));
drawnow;
end

function default_output_3d(par, x, U, step)
clf, imagesc(x{1}, x{2}, U(:, :, ceil(length(x{3})/2))'), axis xy equal tight;
title(sprintf('%s with DG-%d-%d at t = %0.2f', ...
    par.name, par.ord_t, par.ord_x, par.t_plot(step)));
drawnow;
end