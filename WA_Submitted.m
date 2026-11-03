function u = WA_Submitted(EeeT)
% ============================================================
% WA: Error-covariance-weighted average
%
% INPUT
%   EeeT : error covariance matrix (p x p)
%
% OUTPUT
%   u    : merging weight (p x 1)
%
% METHOD
%   The weighted-average solution is obtained by minimizing the variance
%   of the merged estimate under the unbiasedness constraint:
%
%       min_u   u' EeeT u
%       s.t.    eta' u = 1
%
%   where eta is a p x 1 vector of ones.
%
%   The closed-form solution is:
%
%       u = EeeT^(-1) eta / (eta' EeeT^(-1) eta)
%
%   In implementation, MATLAB backslash is used instead of an explicit
%   matrix inverse for numerical stability. If EeeT is nearly singular,
%   a small scale-aware ridge term is added to stabilize the solution.
%
% NOTES
%   - This is the classical variance-minimizing weighted average.
%   - The method uses only the error covariance matrix and does not
%     explicitly account for scaling factors.
%   - For this reason, WA is most directly comparable to other methods
%     when the input datasets are on a common scale.
%
% REFERENCES
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).
%     Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%     IEEE Transactions on Geoscience and Remote Sensing, 60, 1–15.
%
%   Pham, H.T. & Kim, S. (2026).
%     N-Tuple Collocation for Reference-Free Uncertainty Estimation from
%     Multi-Source Hydrological Data: Beyond Triple Collocation and
%     Three-Cornered Hat.
%     Journal of Hydrology (accepted).
%
% If you use the methods presented in the paper and/or this example, 
% please cite this paper where appropriate.
% ============================================================

%% ---------- Input check ----------
if nargin < 1
    error('WA requires one input: EeeT.');
end

[p1,p2] = size(EeeT);

if p1 ~= p2
    error('EeeT must be a square matrix.');
end

%% ---------- Setup ----------
tol = 1e-12;
rcond_min = 1e-10;

% Symmetrize EeeT to suppress tiny numerical asymmetry.
EeeT = 0.5 * (EeeT + EeeT');

if any(isnan(EeeT(:))) || any(isinf(EeeT(:)))
    error('EeeT contains NaN or Inf.');
end

% Remove tiny negative diagonal values caused by round-off.
d = diag(EeeT);
d(d < 0 & d > -tol) = 0;
EeeT(1:p1+1:end) = d;

if any(diag(EeeT) < -tol)
    error('EeeT has negative diagonal entries.');
end

eta = ones(p1,1);

%% ---------- Stabilize near-singular covariance ----------
% Ridge is used only when the covariance matrix is poorly conditioned.
if rcond(EeeT) < rcond_min
    
    scale = mean(diag(EeeT));
    
    if scale <= 0 || isnan(scale) || isinf(scale)
        scale = norm(EeeT,'fro') / max(p1,1);
    end
    
    ridge = max(1e-10 * scale, tol);
    EeeT = EeeT + ridge * eye(p1);
    
end

%% ---------- Solve optimal weights ----------
% Compute EeeT^{-1} eta using backslash.
tmp = EeeT \ eta;

den = eta' * tmp;

% Fallback for extremely ill-conditioned cases.
if abs(den) < tol || isnan(den) || isinf(den)
    tmp = pinv(EeeT) * eta;
    den = eta' * tmp;
end

if abs(den) < tol || isnan(den) || isinf(den)
    error('WA failed: normalization denominator is numerically invalid.');
end

% Normalize so that sum(u) = 1.
u = tmp / den;

% Clean tiny round-off values.
u(abs(u) < tol) = 0;

end
