function data = setup_sldg(data, dim, order, basis_t, poly_t, v, ht, hx)
%SETUP_SLDG Set static data for the SLDG method.
%   This function prepares static data structures required for the
%   Semi-Lagrangian Discontinuous Galerkin (SL DG) method, including
%   including quadrature points, weights, pointwise values and local
%   reference matrices.
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
%      v        - Velocity field.
%      ht       - Time step sizes.
%      hx       - Spatial grid spacing (per dimension).
%
%   Outputs:
%      data - Struct containing precomputed static data for the SLDG method.
%
%   Author: Yi Cai
%   Date: 2024-12-6
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

if dim ~= size(v, 1)
    error('Spatial and velocity dimension dismatches.');
end

%========================================================================
% Initialize parameters
%========================================================================
bbox = [ones(1, dim) * (-1 / 2); ones(1, dim) * (1 / 2)];
deg = order - 1;

if strcmp(poly_t, 'P')
    n_dofs = nchoosek(deg+dim, dim);
elseif strcmp(poly_t, 'Q')
    n_dofs = (deg + 1)^dim;
end

if basis_t == 1 || basis_t == 2
    quad_t = 1;
elseif basis_t == 3
    quad_t = 2;
end

nv = size(v, ndims(v));
nt = numel(ht);
ns = nt * nv;

%========================================================================
% Point of intersection
%========================================================================
r = cell(1, ns);
dt = cell(1, ns);
dx = cell(1, ns);
poi = cell(1, ns);
mtv = multi_index([nt, nv]);
for j = 1:ns
    dt{j} = ht(mtv(j, 1));
    dx{j} = v(:, mtv(j, 2)) .* dt{j};
    r{j} = dx{j} ./ reshape(hx, [], 1);
    poi{j} = r{j} - floor(r{j}) - 1 / 2;
end

%========================================================================
% Volume information
%========================================================================
% Geometry
xc_volp = cell(1, ns);
hx_volp = cell(1, ns);
for j = 1:ns
    xc_volp{j} = zeros(dim, 2^dim);
    hx_volp{j} = zeros(dim, 2^dim);
end
n_volps = zeros(1, ns);
m = linf_index(dim, 2);
for j = 1:ns
    for i = 1:2^dim
        x1 = poi{j};
        x2 = bbox((m(i, :) + (0:dim - 1) * 2)');
        xc_ = (x1 + x2) / 2;
        hx_ = abs(x2-x1);
        if prod(hx_) < 1e-8, continue; end
        n_volps(j) = n_volps(j) + 1;
        xc_volp{j}(:, n_volps(j)) = xc_;
        hx_volp{j}(:, n_volps(j)) = hx_;
    end
    if n_volps(j) == 0
        xc_volp{j} = zeros(dim, 1);
        hx_volp{j} = ones(dim, 1);
        n_volps(j) = 1;
    else
        xc_volp{j} = xc_volp{j}(:, 1:n_volps(j));
        hx_volp{j} = hx_volp{j}(:, 1:n_volps(j));
    end
end

% Quadrature
nq_vol = (deg + 1)^dim;
xq_volp_e = cell(1, ns);
wq_volp_e = cell(1, ns);
xq_volp_u = cell(1, ns);
wq_volp_u = cell(1, ns);
for j = 1:ns
    xq_volp_e{j} = zeros(dim, nq_vol, n_volps(j));
    wq_volp_e{j} = zeros(nq_vol, n_volps(j));
    xq_volp_u{j} = zeros(dim, nq_vol, n_volps(j));
    wq_volp_u{j} = zeros(nq_vol, n_volps(j));
end
for j = 1:ns
    for i = 1:n_volps(j)
        a_e = xc_volp{j}(:, i) - hx_volp{j}(:, i) / 2;
        b_e = xc_volp{j}(:, i) + hx_volp{j}(:, i) / 2;
        a_u = -b_e;
        b_u = -a_e;
        [xq_volp_e{j}(:, :, i), wq_volp_e{j}(:, i)] = ...
            multi_quad(deg+1, a_e, b_e, quad_t, false);
        [xq_volp_u{j}(:, :, i), wq_volp_u{j}(:, i)] = ...
            multi_quad(deg+1, a_u, b_u, quad_t, false);
    end
end

% Basis values
vq_volp_e = cell(1, ns);
vq_volp_u = cell(1, ns);
dvq_volp_e = cell(dim, ns);
for j = 1:ns
    vq_volp_e{j} = zeros(n_dofs, nq_vol, n_volps(j));
    vq_volp_u{j} = zeros(n_dofs, nq_vol, n_volps(j));
    for d = 1:dim
        dvq_volp_e{d, j} = zeros(n_dofs, nq_vol, n_volps(j));
    end
end
for j = 1:ns
    for i = 1:n_volps(j)
        vq_volp_e{j}(:, :, i) = multi_basis( ...
            squeeze(xq_volp_e{j}(:, :, i)), deg, basis_t, poly_t);
        vq_volp_u{j}(:, :, i) = multi_basis( ...
            squeeze(xq_volp_u{j}(:, :, i)), deg, basis_t, poly_t);
        for d = 1:dim
            dvq_volp_e{d, j}(:, :, i) = multi_basis( ...
                squeeze(xq_volp_e{j}(:, :, i)), deg, 1, d, basis_t, poly_t);
        end
    end
end

% Local matrices
v_u_volp = cell(1, ns);
dv_u_volp = cell(dim, ns);
for j = 1:ns
    v_u_volp{j} = zeros(n_dofs, n_dofs, n_volps(j));
    for d = 1:dim
        dv_u_volp{d, j} = zeros(n_dofs, n_dofs, n_volps(j));
    end
end
for j = 1:ns
    for i = 1:n_volps(j)
        % mass
        W = squeeze(wq_volp_e{j}(:, i));
        V = squeeze(vq_volp_e{j}(:, :, i));
        U = squeeze(vq_volp_u{j}(:, :, i));
        v_u_volp{j}(:, :, i) = V * diag(W) * U';
        % gradient
        for d = 1:dim
            V = squeeze(dvq_volp_e{d, j}(:, :, i));
            U = squeeze(vq_volp_u{j}(:, :, i));
            dv_u_volp{j}(:, :, i) = V * diag(W) * U';
        end
    end
end

% Index shift
volp_shift = cell(1, ns);
for j = 1:ns
    volp_shift{j} = zeros(dim, n_volps(j));
end
for j = 1:ns
    mm = ceil(r{j}) .* (r{j} > 0) + (floor(r{j}) + 1) .* (r{j} < 0);
    mm0 = -mm;
    mm1 = -mm + 1;
    for i = 1:n_volps(j)
        bb = xc_volp{j}(:, i) <= 0;
        volp_shift{j}(:, i) = mm0 .* bb + mm1 .* (1 - bb);
    end
end

%========================================================================
% Flux information
%========================================================================
% Geometry
n_flxs = 2 * dim;
n_flxps = zeros(n_flxs, ns);
xc_flxp = cell(n_flxs, ns);
hx_flxp = cell(n_flxs, ns);
for k = 1:ns
    for j = 1:n_flxs
        xc_flxp{j, k} = zeros(dim, 2^(dim - 1));
        hx_flxp{j, k} = zeros(dim, 2^(dim - 1));
    end
end
for k = 1:ns
    for j = 1:n_flxs
        d = floor((j - 1)/2) + 1;
        s = mod(j-1, 2) + 1;
        for i = 1:2^(dim - 1)
            ms = ones(dim, 1) * 2;
            ms(d) = 1;
            m = i2m(i, ms);
            x1 = poi{k};
            x2 = bbox((m + (0:dim - 1) * 2)');
            xc_ = (x1 + x2) / 2;
            hx_ = abs(x2-x1);
            if prod(hx_([1:d-1, d+1:end]))< 1e-8, continue; end
            hx_(d) = 0;
            xc_(d) = bbox(s, d);
            n_flxps(j, k) = n_flxps(j, k) + 1;
            xc_flxp{j, k}(:, n_flxps(j, k)) = xc_;
            hx_flxp{j, k}(:, n_flxps(j, k)) = hx_;
        end
        if n_flxps(j, k) == 0
            xc_flxp{j, k} = zeros(dim, 1);
            xc_flxp{j, k}(d) = bbox(s, d);
            hx_flxp{j, k} = ones(dim, 1);
            hx_flxp{j, k}(d) = 0;
            n_flxps(j, k) = 1;
        else
            xc_flxp{j, k} = xc_flxp{j, k}(:, 1:n_flxps(j, k));
            hx_flxp{j, k} = hx_flxp{j, k}(:, 1:n_flxps(j, k));
        end
    end
end

% Quadrature
nq_flx = (deg + 1)^(dim - 1);
xq_flxp_e = cell(n_flxs, ns);
wq_flxp_e = cell(n_flxs, ns);
xq_flxp_u = cell(n_flxs, ns);
wq_flxp_u = cell(n_flxs, ns);
for k = 1:ns
    for j = 1:n_flxs
        xq_flxp_e{j, k} = zeros(dim, nq_flx, n_flxps(j, k));
        wq_flxp_e{j, k} = zeros(nq_flx, n_flxps(j, k));
        xq_flxp_u{j, k} = zeros(dim, nq_flx, n_flxps(j, k));
        wq_flxp_u{j, k} = zeros(nq_flx, n_flxps(j, k));
    end
end
for k = 1:ns
    for j = 1:n_flxs
        for i = 1:n_flxps(j, k)
            a_e = xc_flxp{j, k}(:, i) - hx_flxp{j, k}(:, i) / 2;
            b_e = xc_flxp{j, k}(:, i) + hx_flxp{j, k}(:, i) / 2;
            a_u = -b_e;
            b_u = -a_e;
            [xq_flxp_e{j, k}(:, :, i), wq_flxp_e{j, k}(:, i)] = ...
                multi_quad(deg+1, a_e, b_e, quad_t, false);
            [xq_flxp_u{j, k}(:, :, i), wq_flxp_u{j, k}(:, i)] = ...
                multi_quad(deg+1, a_u, b_u, quad_t, false);
        end
    end
end

% Basis values
vq_flxp_e_i = cell(n_flxs, ns);
vq_flxp_e_o = cell(n_flxs, ns);
vq_flxp_u_i = cell(n_flxs, ns);
vq_flxp_u_o = cell(n_flxs, ns);
for k = 1:ns
    for j = 1:n_flxs
        vq_flxp_e_i{j, k} = zeros(n_dofs, nq_flx, n_flxps(j, k));
        vq_flxp_e_o{j, k} = zeros(n_dofs, nq_flx, n_flxps(j, k));
        vq_flxp_u_i{j, k} = zeros(n_dofs, nq_flx, n_flxps(j, k));
        vq_flxp_u_o{j, k} = zeros(n_dofs, nq_flx, n_flxps(j, k));
    end
end
for k = 1:ns
    for j = 1:n_flxs
        d = floor((j - 1)/2) + 1;
        s = mod(j-1, 2) + 1;
        for i = 1:n_flxps(j, k)
            ip = repmat({':'}, 1, 3);
            ip{1} = d;
            p = squeeze(xq_flxp_e{j, k}(:, :, i));
            p(ip{:}) = bbox(s, d);
            vq_flxp_e_i{j, k}(:, :, i) = multi_basis(p, deg, basis_t, poly_t);
            p(ip{:}) = bbox(3-s, d);
            vq_flxp_e_o{j, k}(:, :, i) = multi_basis(p, deg, basis_t, poly_t);
            p = squeeze(xq_flxp_u{j, k}(:, :, i));
            p(ip{:}) = -poi{k}(d);
            vq_flxp_u_i{j, k}(:, :, i) = multi_basis(p, deg, basis_t, poly_t);
            if -poi{k}(d) == bbox(1, d)
                p(ip{:}) = bbox(2, d);
                vq_flxp_u_o{j, k}(:, :, i) = multi_basis(p, deg, basis_t, poly_t);
            elseif (-poi{k}(d) == bbox(2, d))
                p(ip{:}) = bbox(1, d);
                vq_flxp_u_o{j, k}(:, :, i) = multi_basis(p, deg, basis_t, poly_t);
            else
                vq_flxp_u_o{j, k}(:, :, i) = multi_basis(p, deg, basis_t, poly_t);
            end
        end
    end
end

% Reference matrices
dv_u_flxp_i = cell(n_flxs, ns);
dv_u_flxp_o = cell(n_flxs, ns);
for k = 1:ns
    for j = 1:n_flxs
        dv_u_flxp_i{j, k} = zeros(n_dofs, n_dofs, n_flxps(j, k));
        dv_u_flxp_o{j, k} = zeros(n_dofs, n_dofs, n_flxps(j, k));
    end
end
for k = 1:ns
    for j = 1:n_flxs
        % gradient
        for i = 1:n_flxps(j, k)
            W = squeeze(wq_flxp_e{j, k}(:, i));
            V = squeeze(vq_flxp_e_i{j, k}(:, :, i));
            U = squeeze(vq_flxp_u_i{j, k}(:, :, i));
            dv_u_flxp_i{j, k}(:, :, i) = V * diag(W) * U';
            U = squeeze(vq_flxp_u_o{j, k}(:, :, i));
            dv_u_flxp_o{j, k}(:, :, i) = V * diag(W) * U';
        end
    end
end

% Index shift
flxp_shift = cell(n_flxs, ns);
for k = 1:nv
    for j = 1:n_flxs
        flxp_shift{j, k} = zeros(dim, n_flxps(j, k));
    end
end
for k = 1:ns
    for j = 1:n_flxs
        for i = 1:n_flxps(j, k)
            mm = ceil(r{k}) .* (r{k} > 0) + (floor(r{k}) + 1) .* (r{k} < 0);
            mm0 = -mm;
            mm1 = -mm + 1;
            bb = xc_flxp{j, k}(:, i) <= 0;
            flxp_shift{j, k}(:, i) = mm0 .* bb + mm1 .* (1 - bb);
        end
    end
end

% Save data
data.nt = nt;
data.nv = nv;
data.ns = ns;
data.poi = poi;
data.dt = dt;
data.dx = dx;

data.n_volps = n_volps;
data.xc_volp = xc_volp;
data.hx_volp = hx_volp;
data.xq_volp_e = xq_volp_e;
data.wq_volp_e = wq_volp_e;
data.xq_volp_u = xq_volp_u;
data.wq_volp_u = wq_volp_u;
data.vq_volp_e = vq_volp_e;
data.vq_volp_u = vq_volp_u;
data.dvq_volp_e = dvq_volp_e;
data.v_u_volp = v_u_volp;
data.dv_u_volp = dv_u_volp;
data.volp_shift = volp_shift;

data.n_flxs = n_flxs;
data.n_flxps = n_flxps;
data.xc_flxp = xc_flxp;
data.hx_flxp = hx_flxp;
data.xq_flxp_e = xq_flxp_e;
data.wq_flxp_e = wq_flxp_e;
data.xq_flxp_u = xq_flxp_u;
data.wq_flxp_u = wq_flxp_u;
data.vq_flxp_e_i = vq_flxp_e_i;
data.vq_flxp_e_o = vq_flxp_e_o;
data.vq_flxp_u_i = vq_flxp_u_i;
data.vq_flxp_u_o = vq_flxp_u_o;
data.dv_u_flxp_i = dv_u_flxp_i;
data.dv_u_flxp_o = dv_u_flxp_o;
data.flxp_shift = flxp_shift;

end
