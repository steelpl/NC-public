function [EeeT_est, theta_est, rho2_est] = TC(covx)
% ============================================================
% Triple Collocation (TC)
% Closed-form estimation assuming mutually uncorrelated errors
%
% INPUT
%   covx : 3x3 covariance matrix of three datasets
%
% OUTPUT
%   EeeT_est  : 3x3 diagonal matrix of estimated error variances
%   theta_est : 3x1 relative scaling vector (theta(1) = 1)
%   rho2_est  : 3x1 squared correlation between each dataset and the latent signal
%
% METHOD SUMMARY
%   Observations: x_i = theta_i * y + e_i,  i = 1,2,3
%   Assumptions:
%     - y is the common true signal with Var(y) = sigma_y^2.
%     - e_i are mutually uncorrelated errors with Var(e_i) = Eii.
%   Population covariance:
%     Q = cov(x) = sigma_y^2 * (theta * theta.') + diag(E11,E22,E33)
%
% ESTIMATION (closed-form for p = 3)
%   Using off-diagonal elements Q12, Q13, Q23 we can eliminate sigma_y^2:
%     E1 = Q11 - (Q12 * Q13) / Q23
%     E2 = Q22 - (Q12 * Q23) / Q13
%     E3 = Q33 - (Q13 * Q23) / Q12
%   These provide estimates of error variances under the assumption that
%   error cross-covariances are zero. A positive estimated sigma_y^2 is
%   obtained from any cyclic relation; we use:
%     sigma_y^2 = (Q12 * Q13) / Q23
%   Relative scaling is fixed by setting theta(1) = 1:
%     theta = [1; Q23/Q13; Q23/Q12]
%   Squared correlations (rho^2) between each x_i and y are:
%     rho2_i = 1 - Eii / sigma_y^2
%
% NOTES
%   - Function only supports p = 3 and requires nonzero off-diagonal
%     covariances (Q12,Q13,Q23). If any are (nearly) zero, returns NaNs.
%   - If any estimated error variance Eii is negative or non-finite the
%     closed-form solution is invalid; returns NaNs.
%   - rho2 is numerically clipped to [0,1] for tiny round-off violations.
%   - This implementation matches the standard Triple Collocation (TC)
%     approach (see Stoffelen, 1998) and is intended to be compact and
%     robust for p = 3.
%
% REFERENCES
%   Stoffelen, A. (1998). Toward the true near‐surface wind speed:
%     Error modeling and calibration using triple collocation.
%     Journal of Geophysical Research: Oceans, 103(C4), 7755-7766.
%   McColl, K. A., et al. (2014). Extended triple collocation:
%     Estimating errors and correlation coefficients with respect to an
%     unknown target. Geophysical Research Letters, 41(17), 6229-6236.


%% ---------- Input check ----------
if nargin < 1
    error('TC requires one input argument: covx');
end
if ~ismatrix(covx) || size(covx,1) ~= size(covx,2)
    error('Input covx must be a square matrix.');
end

p = size(covx,1);
if p ~= 3
    error('TC is only valid for p = 3.');
end

%% ---------- Setup ----------
% Symmetrize to remove tiny asymmetry and extract entries
Q = (covx + covx.') / 2;

Q11 = Q(1,1); Q22 = Q(2,2); Q33 = Q(3,3);
Q12 = Q(1,2); Q13 = Q(1,3); Q23 = Q(2,3);

tol = 1e-12;

%% ---------- Validity check for denominators ----------
% TC requires nonzero cross-covariances
if abs(Q12) < tol || abs(Q13) < tol || abs(Q23) < tol
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

%% ---------- Error variance estimation ----------
% Closed-form TC estimates assuming zero error cross-covariances
E1 = Q11 - (Q12 * Q13) / Q23;
E2 = Q22 - (Q12 * Q23) / Q13;
E3 = Q33 - (Q13 * Q23) / Q12;
Evec = [E1; E2; E3];

%% ---------- Strict validity condition ----------
% Any negative or non-finite error variance invalidates the solution
if any(Evec < 0) || any(~isfinite(Evec))
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

% Diagonal error covariance matrix
EeeT_est = diag(Evec);

%% ---------- Signal variance estimation ----------
% Use cyclic relation; choose the (1,2)-(1,3)/(2,3) formulation
Ey2_est = (Q12 * Q13) / Q23;

if ~(isfinite(Ey2_est) && Ey2_est > 0)
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

%% ---------- Scaling factor estimation ----------
% Normalize to dataset 1 (theta_1 = 1)
theta_est = [1;
             Q23 / Q13;
             Q23 / Q12];

%% ---------- Squared correlation estimation ----------
rho2_est = 1 - Evec / Ey2_est;

% Clip tiny numerical violations to [0,1]
rho2_est(rho2_est < 0 & rho2_est > -tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + tol) = 1;
