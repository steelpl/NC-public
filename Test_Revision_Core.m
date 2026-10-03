function results = Test_Revision_Core(outputDir)
% Focused numerical regression checks for the revised NC/WA primary path.
% Uses base MATLAB only. Run from repository root:
%   Test_Revision_Core(fullfile(pwd,'results','core_checks'));
if nargin < 1
    outputDir = fullfile(fileparts(mfilename('fullpath')),'core_results');
end
if ~exist(outputDir,'dir'); mkdir(outputDir); end
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
rng(83955,'twister');
caseName = strings(0,1); passed = false(0,1); metric = zeros(0,1);

% 1. Default estimator retains decomposition and PSD for full/rank-deficient Q.
worstPSD = 0; worstIdentity = 0; correctionCount = 0;
for k = 1:24
    p = 3 + mod(k,6);
    rankQ = p;
    if mod(k,3) == 0; rankQ = max(1,p-2); end
    A = randn(p,rankQ);
    Q = A*A';
    [E,t,~,ni] = NC(Q,struct('maxIters',400));
    scale = norm(Q,'fro');
    worstPSD = min(worstPSD,min(eig(E))/scale);
    worstIdentity = max(worstIdentity,norm(Q-t*t'-E,'fro')/scale);
    correctionCount = correctionCount + ni.psdRepairApplied;
    assert(ni.feasiblePSD && all(isfinite(t)) && ni.decompositionRelativeError < 1e-13);
end
record("PSD residual across 24 full/rank-deficient covariances",worstPSD > -1e-12,worstPSD);
record("Exact covariance decomposition across 24 cases",worstIdentity < 1e-12,worstIdentity);
record("PSD correction exercised",correctionCount > 0,correctionCount);

% 2. Very large/small covariance units must not change dimensionless results.
Q = [1 .4 .2; .4 1 .3; .2 .3 1];
[E0,t0,r0] = NC(Q);
[w0,~] = WA(E0);
scaleError = 0; weightScaleError = 0;
for scale = [1e-20,1e-8,1e8,1e20]
    [E,t,r] = NC(scale*Q);
    [w,~] = WA(E);
    scaleError = max([scaleError,norm(E/scale-E0,'fro'), ...
        norm(t*t'/scale-t0*t0','fro'),norm(r-r0)]);
    weightScaleError = max(weightScaleError,norm(w-w0));
end
record("NC invariant to covariance unit changes 1e-20 through 1e20",scaleError < 1e-8,scaleError);
record("WA invariant to covariance unit changes 1e-20 through 1e20",weightScaleError < 1e-8,weightScaleError);

% 3. Exact symmetric L1 derivative: infer one solver step independently via
% a central-difference gradient, at a differentiable interior initialization.
[V,D] = eig(Q-.5*eye(3));
[lambda,ix] = max(diag(D));
theta0 = sqrt(lambda)*V(:,ix);
mask = ~eye(3);
loss = @(t) sum(abs(reshape((t*t'-Q).*mask,[],1)));
gradientFD = zeros(3,1); h = 1e-6;
for j = 1:3
    e = zeros(3,1); e(j) = h;
    gradientFD(j) = (loss(theta0+e)-loss(theta0-e))/(2*h);
end
eta0 = 1e-5;
expectedTheta = theta0 - eta0/norm(theta0)*gradientFD;
[~,actualTheta] = NC(Q,struct('enforcePSD',false,'nStarts',1, ...
    'maxIters',1,'minIters',0,'patience',1,'eta0',eta0));
gradientStepError = min(norm(actualTheta-expectedTheta),norm(actualTheta+expectedTheta));
record("One solver step agrees with finite-difference symmetric gradient",gradientStepError < 1e-10,gradientStepError);

% 4. An excessive requested shift cannot produce a spurious zero start.
[~,t,~,ni] = NC(Q,struct('betaFactor',1e8));
record("Adaptive spectral shift retains a nonzero initialization", ...
    ni.spectralShiftAdapted && norm(t)>0 && ni.actualSpectralShift <= .95*max(eig(Q))*(1+1e-14),norm(t));

% 5. Spectral flooring repairs indefinite and singular matrices before solve.
inputs = {[1 2;2 1],diag([1,1e-16,2]),ones(3),[-1 .1;.1 2]};
worstCondition = 0; worstStationarity = 0;
for k = 1:numel(inputs)
    [w,wi] = WA(inputs{k});
    C = wi.covarianceUsed;
    [~,flag] = chol(C);
    assert(flag == 0 && wi.effectiveSPD && wi.regularized && abs(sum(w)-1)<1e-12);
    assert(all(isfinite(w)) && wi.predictedVariance >= 0);
    % First-order optimality of the equality-constrained convex quadratic:
    % C*w must be constant over all components.
    v = C*w;
    stat = norm(v-mean(v)) / max(norm(C,'fro')*norm(w),realmin);
    worstStationarity = max(worstStationarity,stat);
    worstCondition = max(worstCondition,wi.effectiveConditionNumber);
end
record("SPD floor condition bound on four difficult covariances",worstCondition <= 1e8*(1+1e-12),worstCondition);
record("Regularized weights satisfy constrained quadratic optimality",worstStationarity < 1e-12,worstStationarity);

% 6. No unnecessary change for an ordinary SPD covariance.
C = [2 .2;.2 1];
[w,wi] = WA(C);
wExpected = (C\ones(2,1)); wExpected = wExpected/sum(wExpected);
record("Unchanged well-conditioned SPD weighted-average solution", ...
    ~wi.regularized && norm(w-wExpected)<1e-12,norm(w-wExpected));
try
    WA([1 2;2 1],struct('regularize',false));
    rejected = false;
catch ME
    rejected = strcmp(ME.identifier,'WA:NotPositiveDefinite');
end
record("Unregularized solve rejects indefinite covariance",rejected,double(rejected));

% 7. p=2 ambiguity exists even with zero loss and PSD residuals.
Q2 = [1 .3;.3 1];
thetaA = [sqrt(.3);sqrt(.3)]; thetaB = [1;.3];
EA = Q2-thetaA*thetaA'; EB = Q2-thetaB*thetaB';
wA = WA(EA); wB = WA(EB);
[~,~,~,ni] = NC(Q2);
record("p=2 flagged and equally fitting PSD decompositions give different weights", ...
    ni.p2Underdetermined && ~ni.physicalErrorIdentified && ...
    min(eig(EA))>=-1e-14 && min(eig(EB))>=-1e-14 && ...
    abs(EA(1,2))+abs(EB(1,2))<1e-14 && norm(wA-wB)>.5,norm(wA-wB));

% 8. Zero and singular cases have explicit semantics; random starts do not
% overwrite the caller's global random stream.
[EZ,tZ,rZ,nZ] = NC(zeros(3));
[wZ,wzi] = WA(EZ);
record("Zero covariance produces explicit undefined rho2 and equal-weight degeneracy", ...
    all(EZ(:)==0) && all(tZ==0) && all(isnan(rZ)) && nZ.feasiblePSD && ...
    norm(wZ-ones(3,1)/3)==0 && wzi.degenerateZeroCovariance,0);
before = rng;
NC(Q,struct('nRandomStarts',3));
after = rng;
record("Local random starts preserve caller RNG state",isequal(before,after),double(isequal(before,after)));

% 9. Input repair is reported rather than hidden in the decomposition.
[~,~,~,ni] = NC([1 2;2 1]);
record("Indefinite input repair recorded",ni.inputPSDRepaired && ni.inputPSDRelativeChange>0,ni.inputPSDRelativeChange);
try
    NC([1 2;2 1],struct('repairInputPSD',false));
    rejected = false;
catch
    rejected = true;
end
record("Disabling input PSD repair rejects an indefinite input",rejected,double(rejected));

results = table(caseName,passed,metric,'VariableNames',{'Check','Passed','Metric'});
writetable(results,fullfile(outputDir,'core_regression_results.csv'));
save(fullfile(outputDir,'core_regression_results.mat'),'results','worstPSD','worstIdentity', ...
    'correctionCount','scaleError','weightScaleError','gradientStepError', ...
    'worstCondition','worstStationarity','wA','wB');
disp(results);
fprintf('MATLAB %s | %s\n',version,computer);
assert(all(results.Passed),'RevisionCore:Failed','One or more numerical regression checks failed.');

    function record(name,pass,value)
        caseName(end+1,1) = name;
        passed(end+1,1) = pass;
        metric(end+1,1) = value;
    end
end
