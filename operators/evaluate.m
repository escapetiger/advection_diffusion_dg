function res = evaluate(par, ref, u)
%EVALUATE
%   Evaluate DG solution values at prescribed points.

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx);    % Total number of cells
nl = ref.n_dofs;       % Number of local degrees of freedom
np = ref.np;           % Number of points per cell
vp = ref.vp;           % Interpolation matrix (np x nl)
sz = size(u);         % Size of input solution
nb = sz(2:end);       % Additional dimensions beyond the first

%========================================================================
% Vectorized Evaluation
%========================================================================
% Reshape solution to (nl x nc x other dimensions)
u_reshaped = reshape(u, nl, nc, []);

% Compute interpolation for all cells at once
res_reshaped = pagemtimes(vp', u_reshaped);

% Reshape result to match desired output dimensions
res = reshape(res_reshaped, [np * nc, nb]);

end
