function [u, info] = WA(EeeT, opts)
% WA Error-covariance weighted average with an explicit SPD safeguard.
%   u = WA(E) preserves the original one-output calling convention.
%   [u,info] = WA(E,opts) also reports the covariance actually used.
%
%   For a covariance C > 0, minimize u'*C*u subject to sum(u)=1:
%       u = (C\ones(p,1)) / (ones(1,p)*(C\ones(p,1))).
%   Negative weights are permitted; the estimator assumes inputs are on a
%   common unbiased scale. Covariance admissibility does not establish that
%   an estimated covariance equals the physical error covariance.
%
%   The primary path makes E symmetric and replaces every eigenvalue below
%       tau = max(opts.minEigenvalueRatio,p*eps) * max(abs(eig(E)))
%   by tau. The default ratio is 1e-8, bounding spectral conditioning by
%   1e8 (up to floating-point error). This is a spectral floor, not a ridge.
%   Indefinite raw estimates are therefore explicit regularized comparisons.
%   For E=0, all weights minimize the objective; return equal weights and
%   record this degeneracy instead of inventing a unit-dependent floor.
%
%   opts.minEigenvalueRatio : relative eigenvalue floor, default 1e-8
%   opts.regularize        : default true; false requires an SPD input and
%                            is intended for well-conditioned diagnostics
%
%   info.covarianceUsed is in the input covariance units. Important fields:
%     regularized, inputIndefinite, relativeChange, eigenvalueFloor,
%     minEigenvalueBefore, minEigenvalueAfter, effectiveConditionNumber,
%     predictedVariance, effectiveSPD, degenerateZeroCovariance.
%   Changes in weights caused by this safeguard must be reported, not
%   interpreted as evidence that identifiability has been resolved.

if nargin < 1
    error('WA:MissingInput','WA requires an error covariance matrix.');
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isnumeric(EeeT) || ~isreal(EeeT) || ~ismatrix(EeeT) || ...
        isempty(EeeT) || size(EeeT,1) ~= size(EeeT,2) || any(~isfinite(EeeT(:)))
    error('WA:InvalidCovariance','EeeT must be a finite real nonempty square matrix.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error('WA:InvalidOptions','opts must be a scalar structure.');
end
ratio = get_option(opts,'minEigenvalueRatio',1e-8);
regularize = get_option(opts,'regularize',true);
validateattributes(ratio,{'numeric'},{'real','finite','scalar','positive','<=',1},'WA','minEigenvalueRatio');
validateattributes(regularize,{'logical','numeric'},{'scalar','binary'},'WA','regularize');

p = size(EeeT,1);
E = 0.5 * double(EeeT) + 0.5 * double(EeeT');
entryScale = max(abs(E(:)));
info = struct('regularized',false,'inputIndefinite',false, ...
    'relativeChange',0,'eigenvalueFloor',0, ...
    'minEigenvalueBefore',0,'minEigenvalueAfter',0, ...
    'effectiveConditionNumber',Inf,'covarianceUsed',E, ...
    'predictedVariance',0,'effectiveSPD',false, ...
    'degenerateZeroCovariance',entryScale == 0, ...
    'minEigenvalueRatio',ratio,'maxEigenvalueAdjustment',0, ...
    'normalizationError',0,'regularizationMethod','relative spectral floor');
if entryScale == 0
    if ~regularize
        error('WA:NotPositiveDefinite','Unregularized WA requires a strictly positive definite covariance.');
    end
    u = ones(p,1) / p;
    return;
end

% Global normalization makes weights and the safeguard invariant to a
% positive scalar change of covariance units, including very small units.
A = E / entryScale;
[V,D] = eig(A);
lambda = real(diag(D));
lambdaScale = max(abs(lambda));
roundoff = p * eps * lambdaScale;
floorNormalized = max(ratio,p*eps) * lambdaScale;
info.minEigenvalueBefore = min(lambda) * entryScale;
info.inputIndefinite = min(lambda) < -roundoff;
if regularize
    lambdaUsed = max(lambda,floorNormalized);
else
    if min(lambda) <= 0
        error('WA:NotPositiveDefinite','Unregularized WA requires a strictly positive definite covariance.');
    end
    lambdaUsed = lambda;
    floorNormalized = 0;
end
Aused = V * diag(lambdaUsed) * V';
Aused = 0.5 * (Aused + Aused');

% Cholesky establishes positive definiteness before the inverse action.
[R,flag] = chol(Aused);
if flag ~= 0
    error('WA:RegularizationFailure','Covariance is not SPD after the selected safeguard; increase minEigenvalueRatio.');
end
one = ones(p,1);
z = R \ (R' \ one);
denominator = sum(z);
if ~isfinite(denominator) || denominator <= 0 || any(~isfinite(z))
    error('WA:InvalidNormalization','SPD weight normalization failed.');
end
u = z / denominator;
% Preserve the unit-sum constraint; do not threshold individual weights.
u = u / sum(u);

info.regularized = any(lambda < floorNormalized);
info.eigenvalueFloor = floorNormalized * entryScale;
info.maxEigenvalueAdjustment = max(lambdaUsed-lambda) * entryScale;
info.relativeChange = norm(Aused-A,'fro') / norm(A,'fro');
info.covarianceUsed = Aused * entryScale;
info.minEigenvalueAfter = min(eig(Aused)) * entryScale;
info.effectiveConditionNumber = max(lambdaUsed) / min(lambdaUsed);
info.predictedVariance = (u' * Aused * u) * entryScale;
info.effectiveSPD = true;
info.normalizationError = abs(sum(u)-1);
end


function value = get_option(opts,name,defaultValue)
if isfield(opts,name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end
