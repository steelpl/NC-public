function [EeeT_est, theta_est, rho2_est, info] = NC_Submitted(covx, opts)
% ============================================================
% N-Tuple Collocation (NC)
% Constrained rank-1 decomposition of a covariance matrix
%
% INPUT
%   covx : p-by-p covariance matrix
%   opts : optional structure controlling numerical optimization
%          opts.nStarts       : number of deterministic spectral starts
%          opts.nRandomStarts : number of additional random starts
%          opts.randomSeed    : random seed for random starts
%          opts.maxIters      : maximum iterations per start
%          opts.patience      : early-stopping patience
%          opts.enforcePSD    : enforce PSD error covariance, default false
%          opts.repairInputPSD: repair indefinite input covariance, default true
%          opts.verbose       : print optimization summary
%
% OUTPUT
%   EeeT_est  : estimated error covariance (p-by-p), symmetric
%   theta_est : estimated latent-scaling vector (p-by-1)
%   rho2_est  : squared correlation with latent signal (p-by-1)
%   info      : diagnostic structure with objective and convergence details
%
% PURPOSE
%   Find theta and EeeT satisfying
%       Q ≈ theta*theta' + EeeT
%   with theta_i^2 <= Q_ii so that EeeT has nonnegative diagonals. By
%   optionally, the final estimate can be adjusted so that EeeT is positive
%   semidefinite.
%
%   The fit minimizes the off-diagonal L1 mismatch between Q and
%   theta*theta':
%
%       sum_{i ~= j} |Q_ij - theta_i theta_j|
%
% METHOD SUMMARY
%   - Symmetrize input covariance matrix Q.
%   - Initialize theta from several leading eigenpairs of Q - beta*I.
%   - Apply projected subgradient descent with monotone backtracking to
%     reduce the off-diagonal L1 mismatch.
%   - After each step, project theta_i to satisfy |theta_i| <= sqrt(Q_ii).
%   - Stop adaptively when the off-diagonal L1 objective no longer improves.
%   - Return the best theta found across all starts.
%   - If requested, shrink theta only as much as needed so that
%     Q - theta*theta' is positive semidefinite.
%
% NOTES
%   - maxIters is only a safety cap.
%   - Early stopping is controlled by objective improvement, not iteration count.
%   - Small tolerances guard against tiny negative numerical values.
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
%
% If you use the methods presented in the paper and/or this example, 
% please cite this paper where appropriate.
% ============================================================

%% ---------- Input validation ----------
if nargin < 1
    error('NC requires one input argument: covx');
end

if nargin < 2 || isempty(opts)
    opts = struct();
end

if ~ismatrix(covx) || size(covx,1) ~= size(covx,2)
    error('Input covx must be a square matrix.');
end

if any(~isfinite(covx(:)))
    error('Input covx must contain only finite values.');
end

p = size(covx,1);

%% ---------- Parameters ----------
params.betaFactor    = getOpt(opts, 'betaFactor', 0.5);
params.eta0          = getOpt(opts, 'eta0', 0.01);
params.maxIters      = getOpt(opts, 'maxIters', 2000);
params.minIters      = getOpt(opts, 'minIters', 20);
params.patience      = getOpt(opts, 'patience', 30);
params.tol           = getOpt(opts, 'tol', 1e-12);
params.objTol        = getOpt(opts, 'objTol', 1e-10);
params.etaDecay      = getOpt(opts, 'etaDecay', 0.7);
params.etaMin        = getOpt(opts, 'etaMin', 1e-10);
params.lineSearchMax = getOpt(opts, 'lineSearchMax', 12);
params.nStarts       = min(p, max(1, getOpt(opts, 'nStarts', min(4, p))));
params.nRandomStarts = max(0, getOpt(opts, 'nRandomStarts', 0));
params.randomSeed    = getOpt(opts, 'randomSeed', 1);
params.enforcePSD    = getOpt(opts, 'enforcePSD', false);
params.repairInputPSD = getOpt(opts, 'repairInputPSD', true);
params.psdTol        = getOpt(opts, 'psdTol', 1e-10);
params.verbose       = getOpt(opts, 'verbose', false);

%% ---------- Symmetrize and check covariance matrix ----------
Q = 0.5 * (covx + covx');       % remove numerical asymmetry
[Q, inputPSDInfo] = repair_psd_matrix(Q, params.psdTol, params.repairInputPSD, ...
    'Input covariance matrix');
diagQ = diag(Q);

if any(diagQ < -params.tol)
    error('Input covariance matrix has negative diagonal entries.');
end

diagQ(diagQ < 0 & diagQ > -params.tol) = 0;
Q(1:p+1:end) = diagQ;
sqrtDiagQ = sqrt(max(diagQ,0)); % projection bound for theta

beta = params.betaFactor * min(diagQ);  % spectral shift for initialization

%% ---------- Off-diagonal mask ----------
offMask = true(p);
offMask(1:p+1:end) = false;     % exclude diagonal from L1 objective

%% ---------- Multi-start projected subgradient descent ----------
thetaStarts = build_theta_starts(Q,beta,sqrtDiagQ,params);

bestGlobalObj = Inf;
bestGlobalTheta = NaN(p,1);
startInfo = repmat(struct( ...
    'initialObjective', NaN, ...
    'finalObjective', NaN, ...
    'iterations', 0, ...
    'converged', false), size(thetaStarts,2), 1);

for iStart = 1:size(thetaStarts,2)

    [theta_i, obj_i, runInfo] = run_projected_l1( ...
        Q, thetaStarts(:,iStart), sqrtDiagQ, offMask, params);

    startInfo(iStart) = runInfo;

    if obj_i < bestGlobalObj
        bestGlobalObj = obj_i;
        bestGlobalTheta = theta_i;
        bestStart = iStart;
    end
end

%% ---------- Use best solution ----------
theta_est = bestGlobalTheta;

psdInfo = struct('applied', false, 'thetaScale', 1, ...
    'minEigenvalueBefore', NaN, 'minEigenvalueAfter', NaN);
if params.enforcePSD
    [theta_est, psdInfo] = enforce_residual_psd(Q, theta_est, params.psdTol);
end

theta_outer = theta_est * theta_est';
EeeT_est = Q - theta_outer;
EeeT_est = 0.5 * (EeeT_est + EeeT_est');   % enforce symmetry
if params.enforcePSD
    [EeeT_est, finalPSDInfo] = repair_psd_matrix(EeeT_est, params.psdTol, true, ...
        'Estimated error covariance');
else
    finalPSDInfo = psd_matrix_info(EeeT_est);
end

%% ---------- Clean tiny numerical negatives ----------
d = diag(EeeT_est);
d(d < 0 & d > -params.tol) = 0;
EeeT_est(1:p+1:end) = d;

%% ---------- Squared correlation ----------
rho2_est = (theta_est.^2) ./ max(diagQ,params.tol);

rho2_est(rho2_est < 0 & rho2_est > -params.tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + params.tol) = 1;

finalObjective = offdiag_L1_objective(Q, theta_est, offMask);

%% ---------- Diagnostics ----------
info = struct();
info.objective = finalObjective;
info.objectiveBeforePSD = bestGlobalObj;
info.bestStart = bestStart;
info.nStarts = size(thetaStarts,2);
info.startInfo = startInfo;
info.minErrorVariance = min(diag(EeeT_est));
info.minErrorCovEigenvalue = finalPSDInfo.minEigenvalueAfter;
info.maxRho2 = max(rho2_est);
info.feasibleDiagonal = all(diag(EeeT_est) >= -params.tol);
info.feasiblePSD = info.minErrorCovEigenvalue >= -params.psdTol;
info.inputMinEigenvalue = inputPSDInfo.minEigenvalueBefore;
info.inputPSDRepaired = inputPSDInfo.repaired;
info.psdEnforced = params.enforcePSD;
info.psdRepairApplied = psdInfo.applied || finalPSDInfo.repaired;
info.thetaScalePSD = psdInfo.thetaScale;
info.params = params;

if params.verbose
    fprintf('NC: best objective = %.6g, best start = %d/%d, min diag(E) = %.6g, min eig(E) = %.6g\n', ...
        info.objective, info.bestStart, info.nStarts, info.minErrorVariance, ...
        info.minErrorCovEigenvalue);
end

end


%% ============================================================
% Local function 0: Option reader
%% ============================================================
function value = getOpt(opts,name,defaultValue)

if isfield(opts,name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end

end


%% ============================================================
% Local function 1: Initializations
%% ============================================================
function thetaStarts = build_theta_starts(Q,beta,sqrtDiagQ,params)

p = size(Q,1);

[vecs, vals] = eig(Q - beta * eye(p));
eigVals = diag(vals);
[eigVals, order] = sort(eigVals, 'descend');
vecs = vecs(:,order);

nSpectral = min(params.nStarts, p);
thetaStarts = zeros(p, nSpectral + params.nRandomStarts);

for i = 1:nSpectral
    theta = sqrt(max(eigVals(i),0)) * vecs(:,i);
    thetaStarts(:,i) = project_theta(theta,sqrtDiagQ);
end

if params.nRandomStarts > 0
    rng(params.randomSeed);
    for i = 1:params.nRandomStarts
        idx = nSpectral + i;
        theta = sqrtDiagQ .* (2*rand(p,1) - 1);
        thetaStarts(:,idx) = project_theta(theta,sqrtDiagQ);
    end
end

% Remove duplicate zero starts except when every initialization is zero.
norms = vecnorm(thetaStarts);
if any(norms > 0)
    thetaStarts = thetaStarts(:, norms > 0);
end

end


%% ============================================================
% Local function 2: Projected L1 optimization
%% ============================================================
function [bestTheta,bestObj,info] = run_projected_l1(Q,theta0,sqrtDiagQ,offMask,params)

theta = project_theta(theta0,sqrtDiagQ);
obj_prev = offdiag_L1_objective(Q,theta,offMask);

bestTheta = theta;
bestObj = obj_prev;

eta = params.eta0 / max(norm(theta),1e-8);
noImprove = 0;
converged = false;

for k = 1:params.maxIters

    theta_outer = theta * theta';

    % Subgradient of off-diagonal L1 mismatch. The factor of two from the
    % symmetric off-diagonal sum is absorbed into the step-size scale.
    S = sign(theta_outer - Q);
    S(~offMask) = 0;
    grad = S * theta;

    accepted = false;
    etaTry = eta;

    for j = 1:params.lineSearchMax
        theta_new = theta - etaTry * grad;
        theta_new = project_theta(theta_new,sqrtDiagQ);
        obj_new = offdiag_L1_objective(Q,theta_new,offMask);

        if obj_new <= obj_prev || etaTry <= params.etaMin
            accepted = true;
            break;
        end

        etaTry = max(etaTry * params.etaDecay, params.etaMin);
    end

    if ~accepted
        obj_new = obj_prev;
        theta_new = theta;
    end

    eta = etaTry;

    if obj_new < bestObj - params.objTol * max(1,bestObj)
        bestObj = obj_new;
        bestTheta = theta_new;
        noImprove = 0;
    else
        noImprove = noImprove + 1;
    end

    if k >= params.minIters && noImprove >= params.patience
        converged = true;
        break;
    end

    theta = theta_new;
    obj_prev = obj_new;
end

info = struct();
info.initialObjective = offdiag_L1_objective(Q,theta0,offMask);
info.finalObjective = bestObj;
info.iterations = k;
info.converged = converged;

end


%% ============================================================
% Local function 3: Projection
%% ============================================================
function theta = project_theta(theta,sqrtDiagQ)

% Enforce |theta_i| <= sqrt(Q_ii)
theta = sign(theta) .* min(abs(theta),sqrtDiagQ);

end


%% ============================================================
% Local function 4: Off-diagonal L1 objective
%% ============================================================
function obj = offdiag_L1_objective(Q,theta,offMask)

R = theta*theta' - Q;
obj = sum(abs(R(offMask)));

end


%% ============================================================
% Local function 5: PSD handling
%% ============================================================
function [theta, info] = enforce_residual_psd(Q, theta, tol)

info = struct();
info.applied = false;
info.thetaScale = 1;

E0 = Q - theta*theta';
E0 = 0.5 * (E0 + E0');
eigE0 = eig(E0);
info.minEigenvalueBefore = min(eigE0);

scaleQ = max(1, norm(Q, 'fro'));
absTol = tol * scaleQ;
if info.minEigenvalueBefore >= -absTol
    info.minEigenvalueAfter = info.minEigenvalueBefore;
    return;
end

[V, D] = eig(0.5 * (Q + Q'));
qEig = real(diag(D));
qScale = max(1, max(abs(qEig)));
pos = qEig > tol * qScale;
thetaRot = V' * theta;

if any(~pos) && norm(thetaRot(~pos)) > sqrt(tol) * max(norm(theta), 1)
    alpha = 0;
else
    gamma = sum((thetaRot(pos).^2) ./ qEig(pos));
    if gamma <= 1
        alpha = 1;
    else
        alpha = sqrt(1 / gamma);
    end
end

alpha = max(0, min(1, alpha));
theta = alpha * theta;

E1 = Q - theta*theta';
E1 = 0.5 * (E1 + E1');
eigE1 = eig(E1);
info.applied = true;
info.thetaScale = alpha;
info.minEigenvalueAfter = min(eigE1);

end


function [A, info] = repair_psd_matrix(A, tol, allowRepair, labelText)

A = 0.5 * (A + A');
info = psd_matrix_info(A);

scaleA = max(1, norm(A, 'fro'));
absTol = tol * scaleA;
if info.minEigenvalueBefore >= 0
    info.minEigenvalueAfter = info.minEigenvalueBefore;
    return;
end

if ~allowRepair && info.minEigenvalueBefore < -absTol
    error('%s is not positive semidefinite. Minimum eigenvalue: %.6g.', ...
        labelText, info.minEigenvalueBefore);
elseif ~allowRepair
    info.minEigenvalueAfter = info.minEigenvalueBefore;
    return;
end

[V, D] = eig(A);
eigVals = real(diag(D));
eigVals(eigVals < 0) = 0;
A = V * diag(eigVals) * V';
A = 0.5 * (A + A');

minAfter = min(eig(A));
if minAfter < 0
    A = A + (-minAfter + absTol) * eye(size(A));
    A = 0.5 * (A + A');
    minAfter = min(eig(A));
end

info.repaired = true;
info.minEigenvalueAfter = minAfter;

end


function info = psd_matrix_info(A)

eigVals = eig(0.5 * (A + A'));
info = struct();
info.repaired = false;
info.minEigenvalueBefore = min(eigVals);
info.minEigenvalueAfter = info.minEigenvalueBefore;

end
