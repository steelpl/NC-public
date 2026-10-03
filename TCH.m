function [EeeT_est, theta_est, rho2_est, E_TCH_raw] = TCH(covx, ref_idx)
% ============================================================
% Three-Cornered Hat (TCH)
% Reference-free error covariance estimation
%
% INPUT
%   covx    - p-by-p symmetric covariance matrix of observed datasets
%   ref_idx - optional index of reference dataset (integer in 1:p);
%             default is p (last dataset)
%
% OUTPUT
%   EeeT_est  - estimated error covariance matrix (p x p)
%   theta_est - estimated signal standard deviations (p x 1)
%   rho2_est  - squared correlation between each dataset and the latent signal (p x 1)
%   E_TCH_raw - raw TCH error covariance matrix in original ordering (p x p)
%
% NOTES
%   - This function implements a numerically robust Three-Cornered Hat estimator.
%   - It minimizes the quadratic mean of off-diagonal error covariances with
%     penalty terms to enforce positivity and the TCH feasibility constraint.
%   - Uses fminsearch (Nelder–Mead) for unconstrained optimization; options are
%     tuned for high precision. No Optimization Toolbox required.
%   - Returns NaNs for all outputs if the estimated signal variances or error
%     variances are invalid (negative beyond numerical tolerance or non-finite).
%
% REFERENCES
%   Grubbs, F. E. (1948). On estimating precision of measuring instruments
%     and product variability. Journal of the American Statistical Association, 43(242), 243-264.
%   Ferreira, V. G., Montecino, H. D., Yakubu, C. I., & Heck, B. (2016).
%     Uncertainties of GRACE time-variable gravity-field solutions based on the
%     three-cornered hat method. Journal of Applied Remote Sensing, 10(1), 015015.

%% Input validation
% Validate required inputs and set default reference index
if nargin < 1
    error('TCH requires at least one input: covx.');
end
if ~ismatrix(covx) || size(covx,1) ~= size(covx,2)
    error('Input covx must be a square matrix.');
end
p = size(covx,1);
if p < 3
    error('TCH requires at least three datasets.');
end
if nargin < 2 || isempty(ref_idx)
    ref_idx = p;
end
if ~isscalar(ref_idx) || ref_idx < 1 || ref_idx > p || floor(ref_idx) ~= ref_idx
    error('ref_idx must be an integer between 1 and p.');
end

%% Setup and symmetrize
tol = 1e-12;
covx = double(covx);
Q = (covx + covx') / 2;         % Ensure symmetry

% Permute so reference is last
idx_nonref = setdiff(1:p, ref_idx, 'stable');
perm = [idx_nonref, ref_idx];
Qp = Q(perm, perm);
m = p - 1;
u = ones(m,1);

%% Difference covariance S = cov(X_i - X_ref)
S = Qp(1:m,1:m) ...
    - Qp(1:m,p) * ones(1,m) ...
    - ones(m,1) * Qp(p,1:m) ...
    + Qp(p,p) * ones(m,m);
S = (S + S') / 2;  % Symmetrize

% Ensure positive definiteness for inversion
[Ve, De] = eig(S);
d = diag(De);
if any(d <= tol)
    % Replace tiny/nonpositive eigenvalues with tol to stabilize inversion
    d = max(d, tol);
    S = Ve * diag(d) * Ve';
    S = (S + S') / 2;
end

Sinv = inv(S);                % m-by-m inverse
K = det(S)^(1/m);
if ~isfinite(K) || K <= 0
    K = tol;
end

%% Initial guess
sigma0 = zeros(m,1);
varN0  = 1 / (2 * (u' * Sinv * u));   % Heuristic starting noise variance
z0 = double([sigma0; varN0]);

%% Optimization settings
% Use unconstrained search to minimize objective with heavy penalties
options = optimset( ...
    'Display',      'off', ...
    'MaxIter',      5000, ...
    'MaxFunEvals',  20000, ...
    'TolX',         1e-12, ...
    'TolFun',       1e-12);

objfun = @(z) TCH_objective(z, S, Sinv, u, K, tol);
z_est = fminsearch(objfun, z0, options);

sigma = z_est(1:m);
varN  = z_est(m+1);

%% Construct raw TCH error covariance in permuted order
Ehat = S - varN * (u * u') + u * sigma' + sigma * u';
E_raw_perm = [Ehat, sigma; sigma', varN];
E_raw_perm = (E_raw_perm + E_raw_perm') / 2;

% Restore original order
E_TCH_raw = nan(p,p);
E_TCH_raw(perm,perm) = E_raw_perm;

%% Compute extended outputs
Evar   = diag(E_TCH_raw);          % Error variances
theta2 = diag(Q) - Evar;           % theta^2

% Numerical cleanup for tiny negative values
theta2(theta2 < 0 & theta2 > -tol) = 0;

% If invalid estimates, return NaNs
if any(theta2 < 0) || any(Evar < -tol) || any(~isfinite(theta2)) || any(~isfinite(Evar))
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

theta_est = sqrt(max(theta2, 0));

% Error covariance: E = Q - theta*theta' but enforce diagonals to Evar
EeeT_est = Q - theta_est * theta_est';
EeeT_est(1:p+1:end) = Evar;
EeeT_est = (EeeT_est + EeeT_est') / 2;

rho2_est = theta2 ./ diag(Q);
rho2_est(rho2_est < 0 & rho2_est > -tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + tol) = 1;
end

%% Local objective function (optimized and well-formatted)
function f = TCH_objective(z, S, Sinv, u, K, tol)
% Objective for TCH: minimize quadratic mean of off-diagonal error covariances
% with penalties enforcing constraints and positive diagonals.

m = length(u);
sigma = z(1:m);
varN  = z(m+1);

% Invalid parameters -> large cost
if varN <= 0 || any(~isfinite(z))
    f = 1e30;
    return
end

% Reconstruct E (permuted)
Ehat = S - varN * (u * u') + u * sigma' + sigma * u';
E = [Ehat, sigma; sigma', varN];
E = (E + E') / 2;

% Objective: quadratic mean of off-diagonal elements, normalized by K^2
off = E - diag(diag(E));
F = sum(off(:).^2) / max(K^2, tol);

% TCH positive-definiteness constraint:
q = sigma - varN * u;
constraint_value = varN - q' * Sinv * q;

penalty = 0;
if constraint_value <= tol
    % Strong penalty when constraint near or violated
    penalty = penalty + 1e12 * (tol - constraint_value)^2 + 1e8;
end

% Penalize significantly negative diagonals
d = diag(E);
neg = d(d < -tol);
if ~isempty(neg)
    penalty = penalty + 1e12 * sum(neg.^2);
end

f = F + penalty;
end