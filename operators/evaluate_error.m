function res = evaluate_error(par, ref, u, x, hx, t)
%EVALUATE_ERROR
%   Compute the L1, L2, and Linf error norms for a DG solution.

%========================================================================
% Parse parameters
%========================================================================
nc = prod(par.nx);   % Total number of cells
nl = ref.n_dofs;      % Local degrees of freedom
xr = ref.xr;          % Reference points for evaluation
wr = ref.wr;          % Weights for integration
vr = ref.vr;          % Evaluation matrix
h = prod(hx);        % Cell volume

%========================================================================
% Precompute Evaluation Points
%========================================================================
% Compute all multi-dimensional cell indices
m = multi_index(par.nx);

% Compute physical evaluation points for all cells
X = cell(1, par.dim);
for d = 1:par.dim
    % X{d} is (nr x nc): evaluation points along dimension `d`
    X{d} = x{d}(m(:, d)) + xr(d, :)' * hx(d);
end

% Evaluate exact solution for all points
Ue = par.fn_exact(par, X, t); % (nr x nc)

%========================================================================
% Evaluate DG Solution
%========================================================================
% Reshape `u` for all cells
u_reshaped = reshape(u, nl, nc, []);

% Evaluate DG solution at reference points for all cells
U = pagemtimes(vr', u_reshaped); % (nr x nc x nb)

%========================================================================
% Compute Residual Norms
%========================================================================
% Compute errors
R = abs(Ue - U); % (nr x nc x nb)

% Compute L1, L2, and Linf norms
res = struct;
res.L1 = h * sum(wr * R, 'all');
res.L2 = sqrt(h * sum(wr * R.^2, 'all'));
res.Linf = max(R, [], 'all');

end
