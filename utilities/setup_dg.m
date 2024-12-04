function [data] = setup_dg(dim, order, basis_t, poly_t, n_plot, n_error)
%SETUP_DG Set static data for the Discontinuous Galerkin (DG) method.
%   This function prepares static data structures required for the
%   Discontinuous Galerkin (DG) method, including quadrature points,
%   weights, basis functions, and derivative matrices.
%
%   Syntax:
%      data = SETUP_DG(dim, order, basis, poly_t)
%
%   Inputs:
%      dim     - Number of spatial dimensions (positive integer).
%      order   - Polynomial order of approximation (positive integer).
%      basis_t - Basis type:
%                  1 = Legendre basis
%                  2 = Chebyshev basis
%                  3 = Lobatto basis
%      poly_t  - Polynomial type:
%                  'P' for simplex polynomials
%                  'Q' for tensor-product polynomials
%      n_plot  - Number of points to plot
%      n_error  - Number of points used to compute residual
%
%   Outputs:
%      data - A struct containing static data for the DG method.
%
%   Example:
%      data = setup_dg(2, 3, 1, 'P');
%
%   Author: Yi Cai
%   Date: 2024-11-21
%   Version: 1.1

%========================================================================
% Check inputs
%========================================================================
if dim < 1 || order < 1 || ~ismember(basis_t, [1, 2, 3]) || ...
        ~ismember(poly_t, {'P', 'Q'})
    error('Invalid input arguments. Please check documentation.');
end

if basis_t == 3 && order == 1
    error('The order of Gauss-Lobatto nodal DG method must > 1.');
end

if dim > 1 && basis_t > 1 && poly_t == 'P'
    error('Nodal DG method requires tensor-product polynomials.');
end

if nargin < 5
    n_plot = 1;
end

if nargin < 6
    n_error = 0;
end


%========================================================================
% Initialize parameters
%========================================================================
bbox = [ones(1, dim) * (-1/2); ones(1, dim) * (1/2)];
deg = order - 1;

if strcmp(poly_t, 'P')
    n_dofs = nchoosek(deg + dim, dim);
elseif strcmp(poly_t, 'Q')
    n_dofs = (deg + 1)^dim;
end

if basis_t == 1 || basis_t == 2
    quad_t = 1;
elseif basis_t == 3
    quad_t = 2;
end

%========================================================================
% Volume information
%========================================================================
% Quadrature
nq_vol = (deg + 1)^dim;
[xq_vol, wq_vol] = multi_quad(deg + 1, bbox(1, :), bbox(2, :), quad_t, false);

% Basis values at quadrature points
vq_vol = multi_basis(xq_vol, deg, basis_t, poly_t);
dvq_vol = cell(1, dim);
for d = 1:dim
    dvq_vol{d} = multi_basis(xq_vol, deg, 1, d, basis_t, poly_t);
end

% Mass matrices
v_u_vol = vq_vol * diag(wq_vol) * vq_vol';

% Derivative matrices
dv_u_vol = cell(1, dim);
for d = 1:dim
    dv_u_vol{d} = dvq_vol{d} * diag(wq_vol) * vq_vol';
end

%========================================================================
% Flux information
%========================================================================
n_flxs = 2 * dim;

% Quadrature
nq_flx = (deg + 1)^(dim - 1);
xq_flx = cell(1, n_flxs);
wq_flx = cell(1, n_flxs);
for d = 1:dim
    i = 2 * d - 1; 
    a = ones(1, dim) * (-1/2);
    b = ones(1, dim) * (1/2);
    a(d) = -1/2; b(d) = -1/2;
    [xq_flx{i}, wq_flx{i}] = multi_quad(deg + 1, a, b, quad_t);
    
    i = 2 * d;
    a = ones(1, dim) * (-1/2);
    b = ones(1, dim) * (1/2);
    a(d) = 1/2; b(d) = 1/2;
    [xq_flx{i}, wq_flx{i}] = multi_quad(deg + 1, a, b, quad_t);
end

% Basis values at quadrature points
vq_flx = cell(1, n_flxs);
for i = 1:n_flxs
    vq_flx{i} = multi_basis(xq_flx{i}, deg, basis_t, poly_t);
end

% Flux matrices
dv_u_flx_i = cell(1, n_flxs);
dv_u_flx_o = cell(1, n_flxs);
for d = 1:dim
    i1 = 2 * d - 1; i2 = 2 * d;
    W1 = diag(wq_flx{i1});
    W2 = diag(wq_flx{i2});
    V1 = squeeze(vq_flx{i1});
    V2 = squeeze(vq_flx{i2});
    dv_u_flx_i{i1} = V1 * W1 * V1';
    dv_u_flx_o{i1} = V1 * W1 * V2';
    dv_u_flx_i{i2} = V2 * W2 * V2';
    dv_u_flx_o{i2} = V2 * W2 * V1';
end

%========================================================================
% Plotting information
%========================================================================
if n_plot > 0
    np = n_plot^dim;
    xp = cell(1, dim);
    for d = 1:dim
        xp{d} = linspace(-1/2, 1/2, n_plot+2);
        xp{d} = xp{d}(2:end - 1);
    end
    [xp_vol{1:dim}] = ndgrid(xp{:});
    xp_vol = reshape(cat(dim+1, xp_vol{:}), [], dim).';
    vp = multi_basis(xp_vol, deg, basis_t, poly_t);
end

%========================================================================
% Error computation
%========================================================================
if n_error > 0
    nr = n_error^dim;
    [xr, wr] = multi_quad(n_error, bbox(1, :), bbox(2, :), quad_t, false);
    vr = multi_basis(xr, deg, basis_t, poly_t);
end

% Assemble data structure
data = struct;
data.n_dofs = n_dofs;
data.nq_vol = nq_vol;
data.xq_vol = xq_vol;
data.wq_vol = wq_vol;
data.vq_vol = vq_vol;
data.dvq_vol = dvq_vol;
data.v_u_vol = v_u_vol;
data.dv_u_vol = dv_u_vol;
data.n_flxs = n_flxs;
data.nq_flx = nq_flx;
data.xq_flx = xq_flx;
data.wq_flx = wq_flx;
data.vq_flx = vq_flx;
data.dv_u_flx_i = dv_u_flx_i;
data.dv_u_flx_o = dv_u_flx_o;

if n_plot > 0
    data.np = np;
    data.xp = xp;
    data.vp = vp;
end

if n_error > 0
    data.nr = nr;
    data.xr = xr;
    data.wr = wr;
    data.vr = vr;
end


end
