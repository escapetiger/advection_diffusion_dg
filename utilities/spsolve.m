function x = spsolve(A, b, varargin)
% SPSOLVE Selects the best solver for a given matrix and right-hand side.
%
%   x = SPSOLVE(A, b) solves the system Ax = b by selecting the most
%   suitable solver based on the properties of the matrix A.
%
%   x = SPSOLVE(A, b, options) allows additional options to be passed
%   as a struct with fields:
%       - tol: Tolerance for iterative solvers (default: 1e-6)
%       - maxit: Maximum iterations for iterative solvers (default: 500)
%       - size_threshold: Matrix size threshold for iterative solvers
%                        (default: 100,000 non-zeros)
%       - is_spd: Logical flag for symmetric positive definite property
%                 (default: [])
%
%   Inputs:
%       - A: Coefficient matrix (sparse or dense)
%       - b: Right-hand side vector
%       - options: Struct for solver options (optional)
%
%   Output:
%       - x: Solution to the system Ax = b

%=========================================================================
% Parse inputs
%=========================================================================
if nargin < 3
    options = struct();
else
    options = varargin{1};
end

% Set default options
% Tolerance for iterative solvers (default: 1e-6)
tol = 1e-6;
if isfield(options, 'tol'), tol = options.tol; end

% Maximum iterations for iterative solvers (default: 500)
maxit = 500;
if isfield(options, 'maxit'), maxit = options.maxit; end

% Size threshold for selecting iterative solvers (default: 100,000 non-zeros)
size_threshold = 1e5;
if isfield(options, 'size_threshold'), size_threshold = options.size_threshold; end

% User-provided flag for symmetric positive definiteness (default: empty)
is_spd = [];
if isfield(options, 'is_spd'), is_spd = options.is_spd; end

% Check if A is sparse
is_sparse = issparse(A);

%=========================================================================
% Analyze matrix properties
%=========================================================================
% Ensure the matrix is square
is_square = size(A, 1) == size(A, 2);
if ~is_square
    error('Matrix A must be square.');
end

% Check symmetry and positive definiteness only if not provided by the user
if isempty(is_spd)
    if ~issymmetric(A)
        is_spd = false;
    else
        [~, p] = chol(A);
        is_spd = (p == 0);
    end
end

% Determine the size of the matrix
if is_sparse
    nnz_A = nnz(A); % Number of non-zero entries
else
    nnz_A = numel(A); % Total number of entries for dense matrices
end

%=========================================================================
% Select solver based on matrix properties
%=========================================================================
if nnz_A > size_threshold
    % Large matrix: Use iterative solvers
    if is_spd
        solver = 'pcg'; % Preconditioned Conjugate Gradient for SPD matrices
    elseif is_sparse
        solver = 'gmres'; % Generalized Minimum Residual Method for sparse matrices
    else
        solver = 'bicgstab'; % Bi-Conjugate Gradient Stabilized Method for dense matrices
    end
else
    % Small matrix: Use direct solvers
    if is_spd
        solver = 'chol'; % Cholesky factorization for SPD matrices
    elseif is_sparse
        solver = 'backslash'; % Sparse LU/QR solver
    else
        solver = 'backslash'; % Dense LU solver
    end
end

%=========================================================================
% Solve the system using the selected solver
%=========================================================================
switch solver
    case 'chol'
        % Use Cholesky factorization for SPD matrices
        R = chol(A);
        x = R \ (R' \ b);
    case 'pcg'
        % Use Preconditioned Conjugate Gradient for large SPD matrices
        x = pcg(A, b, tol, maxit);
    case 'gmres'
        % Use Generalized Minimum Residual Method for general sparse matrices
        x = gmres(A, b, [], tol, maxit);
    case 'bicgstab'
        % Use Bi-Conjugate Gradient Stabilized Method for general dense matrices
        x = bicgstab(A, b, tol, maxit);
    case 'backslash'
        % Use MATLAB's backslash operator for direct solvers
        x = A \ b;
    otherwise
        error('No suitable solver found.');
end

end
