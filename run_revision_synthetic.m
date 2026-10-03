function out = run_revision_synthetic(nMC, outDir, completedRuntimeSeconds)
%RUN_REVISION_SYNTHETIC Paired JoH revision experiments (base MATLAB only).
% Run from the folder containing this file.
%   out = run_revision_synthetic(100);
% Smoke: run_revision_synthetic(3, fullfile(pwd,'results','smoke','main'));
% All methods see the same sample covariance; truth never enters an estimator.
% p=3 and unit gains permit a fair, directly comparable TC/TCH/NC benchmark.
% Seasonal and stationary cases use IDENTICAL error draws when alpha=0.
% Main merging metrics apply the SAME revised WA SPD stabilization to every
% covariance estimate. UnsafeMergeMSE is a diagnostic of the old unprotected
% solve; it is never the recommended merging result.

if nargin < 1 || isempty(nMC), nMC = 100; end
finalizeOnly=ischar(nMC) && strcmp(nMC,'finalize_existing');
if finalizeOnly, nMC=NaN; end
if nargin<3, completedRuntimeSeconds=NaN; end
here = fileparts(mfilename('fullpath'));
codeDir = here;
if nargin < 2 || isempty(outDir), outDir = fullfile(here,'results'); end
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath(codeDir);
assert(finalizeOnly || (nMC >= 2 && nMC == fix(nMC)), 'nMC must be an integer >= 2.');
started = tic;
cfg = struct('seed',20260917,'nMC',nMC,'p',3,'nTrain',1095,'nTest',1095, ...
    'burnIn',1000,'errorAR',0.5,'stationaryTruthAR',0.8, ...
    'slowTruthAR',0.95,'fastTruthAR',0.4,'annualPeriod',365, ...
    'SNRdB',[0.1 3 10],'baseErrorCorrelation',[0 0.4 0.99], ...
    'signalErrorCorrelation',[0 0.3],'relativeErrorVariance',[0.5 1 1.5], ...
    'truthTypes',{{'stationary_AR','seasonal_multiscale'}}, ...
    'gains',[1 1 1],'ciMethod','Monte Carlo mean +/- 1.96 SE; rates Wilson 95%', ...
    'ncOptions',struct('nStarts',3,'nRandomStarts',0), ...
    'mergeRule','Identical default revised WA PSD projection + SPD floor for every estimator');
methods = {'TC','TCH','NC_raw','NC_PSD','equal_weight','oracle'};
scenarios = make_scenarios(cfg);
nScen = height(scenarios);
if ~finalizeOnly
blank = blank_record();
records = repmat(blank,nMC*nScen*numel(methods),1);
ir = 0;
fprintf('Revision experiment: %d replicates x %d scenarios; n=%d train/%d test.\n',nMC,nScen,cfg.nTrain,cfg.nTest);
for k = 1:nMC
    % Reset once per replicate: scenario order cannot change paired samples.
    rng(cfg.seed+k,'twister');
    train = base_draws(cfg.nTrain,cfg);
    test = base_draws(cfg.nTest,cfg);
    for s = 1:nScen
        sc = scenarios(s,:);
        snr = 10^(sc.SNRdB/10);
        v = cfg.relativeErrorVariance(:)/snr;
        rho = sc.ErrorCorrelation;
        alpha = sc.SignalErrorCorrelation;
        R0 = (1-rho)*eye(cfg.p)+rho*ones(cfg.p);
        L0 = chol(R0,'lower');
        Rtrue = (1-alpha^2)*R0 + alpha^2*ones(cfg.p);
        Etrue = diag(sqrt(v))*Rtrue*diag(sqrt(v));
        idxTruth = sc.TruthIndex;
        ytr = train.truth(:,idxTruth); yte = test.truth(:,idxTruth);
        etr = (sqrt(1-alpha^2)*(train.error*L0') + alpha*ytr*ones(1,cfg.p)).*sqrt(v)';
        ete = (sqrt(1-alpha^2)*(test.error*L0') + alpha*yte*ones(1,cfg.p)).*sqrt(v)';
        Xtr = ytr + etr;
        Xte = yte + ete;
        Q = cov(Xtr);
        Eemp = cov(etr);
        CyE = cov([ytr etr]);
        realizedCorr = CyE(1,2:end)./sqrt(diag(CyE(2:end,2:end)))';
        estimatorSecs = zeros(1,4);
        matrices = cell(1,4); ni = cell(1,4); problem = strings(1,4); tchFit=nan(cfg.p);
        for m = 1:4
            estStart = tic;
            try
                switch m
                    case 1, matrices{m} = TC(Q);
                    case 2, [matrices{m},~,~,tchFit] = TCH(Q);
                    case {3,4}
                        opts = cfg.ncOptions; opts.enforcePSD = (m==4);
                        [matrices{m},~,~,ni{m}] = NC(Q,opts);
                end
            catch err
                matrices{m} = nan(cfg.p); problem(m) = string(err.identifier)+": "+string(err.message);
            end
            estimatorSecs(m) = toc(estStart);
        end
        for m = 1:numel(methods)
            ir = ir+1; r = blank;
            r.Scenario = s; r.Replicate = k; r.Seed = cfg.seed+k;
            r.Truth = string(sc.Truth{1}); r.SNRdB = sc.SNRdB;
            r.ErrorCorrelation = rho; r.SignalErrorCorrelation = alpha;
            r.Method = string(methods{m});
            r.RealizedSignalErrorCorrelation = mean(realizedCorr);
            r.TargetErrorVariance1=v(1); r.TargetErrorVariance2=v(2); r.TargetErrorVariance3=v(3);
            r.TargetMeanErrorVariance=mean(v);
            r.SampleTruthVariance=var(ytr);
            r.OracleMergeMSE = NaN;
            if m <= 4
                E = matrices{m}; r.EstimatorSeconds = estimatorSecs(m); r.EstimatorError = problem(m);
                if m==2 && all(isfinite(tchFit(:)))
                    r.TCHFittedMinEigenvalue=min(eig((tchFit+tchFit')/2));
                    r.TCHFittedNegativeVariance=any(diag(tchFit)<-1e-9*norm(Q,2));
                    r.TCHReconstructionRelativeChange=norm(E-tchFit,'fro')/max(norm(tchFit,'fro'),realmin);
                end
                r.NonfiniteEstimate = any(~isfinite(E(:)));
                if ~r.NonfiniteEstimate
                    E=(E+E')/2; dv = diag(E)-v;
                    r.EstimatedVariance1=E(1,1); r.EstimatedVariance2=E(2,2); r.EstimatedVariance3=E(3,3);
                    r.VarianceBias = mean(dv);
                    r.VarianceMSE = mean(dv.^2);
                    r.NormalizedVarianceMSE = r.VarianceMSE/mean(v)^2;
                    r.EmpiricalVarianceMSE = mean((diag(E)-diag(Eemp)).^2);
                    r.NormalizedCovarianceError = norm(E-Etrue,'fro')/norm(Etrue,'fro');
                    r.MinEigenvalue = min(eig(E));
                    tol=1e-9*max(norm(Q,2),realmin);
                    r.NegativeVariance = any(diag(E)<-tol);
                    r.IndefiniteCovariance = r.MinEigenvalue < -tol;
                    r.ConditionNumber = cond(E);
                    r.UnsafeMergeMSE = unsafe_merge_mse(E,Xte,yte);
                    if ~r.NegativeVariance
                        try
                            [w,wi]=WA(E);
                            r=merge_metrics(r,w,Xte,yte,mean(v));
                            r.MergeCovarianceRelativeChange=get_numeric_field(wi,{'relativeChange','relativeFrobeniusChange','relativeAdjustment','relativeCovarianceChange'});
                            r.MergeFloorApplied=get_numeric_field(wi,{'regularized','floorApplied','spdFloorApplied','regularizationApplied'});
                        catch err
                            r.MergeError=string(err.identifier)+": "+string(err.message);
                        end
                    end
                end
                if ~isempty(ni{m})
                    r.NCObjective=ni{m}.objective;
                    r.NCObjectiveBeforePSD=ni{m}.objectiveBeforePSD;
                    r.NCThetaScale=ni{m}.thetaScalePSD;
                    r.NCPSDRepairApplied=ni{m}.psdRepairApplied;
                    r.NCIterations=sum([ni{m}.startInfo.iterations]);
                end
            elseif m==5
                r=merge_metrics(r,ones(cfg.p,1)/cfg.p,Xte,yte,mean(v));
            else
                [w,~]=WA(Etrue);
                r=merge_metrics(r,w,Xte,yte,mean(v));
            end
            records(ir)=r;
        end
    end
    if mod(k,max(1,round(nMC/10)))==0 || k==1
        fprintf('%d/%d replicates complete; %.1f seconds elapsed.\n',k,nMC,toc(started));
    end
end
replicates=struct2table(records);
summary=summarize_records(replicates,scenarios,methods);
paired=paired_temporal_summary(replicates,scenarios,methods);
comparisons=common_subset_comparisons(replicates,scenarios);
runtimeSeconds=toc(started);
else
    % Recover artifact export after a graphics/file-export interruption.
    % Numerical experiments are not rerun or altered; persisted CSV metrics
    % are loaded, with CSV precision, and this fact is recorded in runInfo.
    replicates=readtable(fullfile(outDir,'replicate_metrics.csv'),'TextType','string');
    % readtable can infer an entirely empty text column as numeric NaNs.
    for exceptionColumn={'EstimatorError','MergeError'}
        label=exceptionColumn{1}; x=string(replicates.(label));
        x(ismissing(x) | x=="NaN")=""; replicates.(label)=x;
    end
    summary=summarize_records(replicates,scenarios,methods);
    paired=readtable(fullfile(outDir,'paired_seasonal_minus_stationary.csv'),'TextType','string');
    comparisons=readtable(fullfile(outDir,'common_subset_method_comparisons.csv'),'TextType','string');
    cfg.nMC=numel(unique(replicates.Replicate)); nMC=cfg.nMC;
    runtimeSeconds=completedRuntimeSeconds;
end
runInfo=struct('matlabVersion',version,'computer',computer,'timestampUTC', ...
    char(datetime('now','TimeZone','UTC')),'runtimeSeconds',runtimeSeconds, ...
    'runtimeDefinition','Generation, fitting and numeric summaries; excludes file/figure export', ...
    'exportRecoveredFromCSV',finalizeOnly, ...
    'NCPath',which('NC'),'WAPath',which('WA'),'TCPath',which('TC'),'TCHPath',which('TCH'), ...
    'toolboxes',ver,'numberOfCovarianceFits',nMC*nScen*4);
writetable(replicates,fullfile(outDir,'replicate_metrics.csv'));
writetable(summary,fullfile(outDir,'scenario_summary.csv'));
writetable(paired,fullfile(outDir,'paired_seasonal_minus_stationary.csv'));
writetable(comparisons,fullfile(outDir,'common_subset_method_comparisons.csv'));
writetable(scenarios,fullfile(outDir,'scenario_configuration.csv'));
sourceCode=struct('NC',fileread(which('NC')),'WA',fileread(which('WA')), ...
    'TC',fileread(which('TC')),'TCH',fileread(which('TCH')),'runner',fileread([mfilename('fullpath') '.m']));
save(fullfile(outDir,'revision_synthetic_results.mat'),'cfg','runInfo','scenarios','replicates','summary','paired','comparisons','sourceCode','-v7');
fid=fopen(fullfile(outDir,'run_configuration.json'),'w');
assert(fid>=0); fprintf(fid,'%s',jsonencode(struct('config',cfg,'runInfo',runInfo),'PrettyPrint',true)); fclose(fid);
% Document figures are generated once by Build_Synthetic_Revision_Artifacts.
write_run_note(cfg,runInfo,summary,paired,outDir);
out=struct('config',cfg,'runInfo',runInfo,'summary',summary,'paired',paired,'outDir',outDir);
fprintf('Completed in %.1f seconds. Results: %s\n',runtimeSeconds,outDir);
end

function scenarios=make_scenarios(cfg)
s=0; rec=struct([]);
for t=1:2
    for snr=cfg.SNRdB
        for rho=cfg.baseErrorCorrelation
            for alpha=cfg.signalErrorCorrelation
                s=s+1;
                rec(s).Scenario=s; rec(s).TruthIndex=t; rec(s).Truth=cfg.truthTypes{t};
                rec(s).SNRdB=snr; rec(s).ErrorCorrelation=rho; rec(s).SignalErrorCorrelation=alpha;
            end
        end
    end
end
scenarios=struct2table(rec);
end

function b=base_draws(n,cfg)
nAll=n+cfg.burnIn;
b.error=ar_draw(randn(nAll,cfg.p),cfg.errorAR,cfg.burnIn);
z=randn(nAll,3);
y1=ar_draw(z(:,1),cfg.stationaryTruthAR,cfg.burnIn);
slow=ar_draw(z(:,2),cfg.slowTruthAR,cfg.burnIn);
fast=ar_draw(z(:,3),cfg.fastTruthAR,cfg.burnIn);
t=(1:n)';
annual=sqrt(2)*sin(2*pi*t/cfg.annualPeriod);
y2=sqrt(0.5)*annual+sqrt(0.3)*slow+sqrt(0.2)*fast;
b.truth=[unit_variance(y1),unit_variance(y2)];
end

function x=ar_draw(z,phi,burn)
x=filter(sqrt(1-phi^2),[1 -phi],z);
x=x(burn+1:end,:);
end

function x=unit_variance(x)
x=(x-mean(x))/std(x,0);
end

function r=blank_record()
r=struct('Scenario',NaN,'Replicate',NaN,'Seed',NaN,'Truth',"",'SNRdB',NaN, ...
    'ErrorCorrelation',NaN,'SignalErrorCorrelation',NaN,'Method',"", ...
    'RealizedSignalErrorCorrelation',NaN,'TargetErrorVariance1',NaN, ...
    'TargetErrorVariance2',NaN,'TargetErrorVariance3',NaN,'TargetMeanErrorVariance',NaN, ...
    'SampleTruthVariance',NaN,'EstimatedVariance1',NaN,'EstimatedVariance2',NaN, ...
    'EstimatedVariance3',NaN,'VarianceBias',NaN,'VarianceMSE',NaN, ...
    'NormalizedVarianceMSE',NaN,'EmpiricalVarianceMSE',NaN,'NormalizedCovarianceError',NaN, ...
    'NonfiniteEstimate',false,'NegativeVariance',false,'IndefiniteCovariance',false, ...
    'MinEigenvalue',NaN,'ConditionNumber',NaN,'MergeMSE',NaN,'NormalizedMergeMSE',NaN, ...
    'MergeR2',NaN,'Weight1',NaN,'Weight2',NaN,'Weight3',NaN,'WeightL1',NaN, ...
    'MaxAbsWeight',NaN,'UnsafeMergeMSE',NaN,'OracleMergeMSE',NaN, ...
    'MergeCovarianceRelativeChange',NaN,'MergeFloorApplied',NaN, ...
    'NCObjective',NaN,'NCObjectiveBeforePSD',NaN,'NCThetaScale',NaN, ...
    'NCPSDRepairApplied',NaN,'NCIterations',NaN,'EstimatorSeconds',NaN, ...
    'TCHFittedMinEigenvalue',NaN,'TCHFittedNegativeVariance',NaN,'TCHReconstructionRelativeChange',NaN, ...
    'EstimatorError',"",'MergeError',"");
end

function r=merge_metrics(r,w,Xte,yte,vmean)
yhat=Xte*w;
r.MergeMSE=mean((yhat-yte).^2);
r.NormalizedMergeMSE=r.MergeMSE/vmean;
C=corrcoef(yhat,yte); r.MergeR2=C(1,2)^2;
r.Weight1=w(1); r.Weight2=w(2); r.Weight3=w(3);
r.WeightL1=sum(abs(w)); r.MaxAbsWeight=max(abs(w));
end

function mse=unsafe_merge_mse(E,X,y)
% Explicit historical finite-weight diagnostic, NOT a covariance-valid rule.
% Skip a singular solve; report failures rather than hide them in a fallback.
mse=NaN;
if any(~isfinite(E(:))) || rcond(E)<1e-12, return; end
a=E\ones(size(E,1),1); den=sum(a);
if any(~isfinite(a)) || ~isfinite(den) || abs(den)<1e-12, return; end
w=a/den;
mse=mean((X*w-y).^2);
end

function v=get_numeric_field(s,names)
v=NaN;
for k=1:numel(names)
    if isfield(s,names{k}), v=double(s.(names{k})); return; end
end
end

function out=summarize_records(T,S,methods)
recs=struct([]); k=0;
for s=1:height(S)
    for m=1:numel(methods)
        D=T(T.Scenario==s & T.Method==methods{m},:);
        k=k+1; r=table2struct(S(s,:)); r.Method=string(methods{m}); r.NReplicates=height(D);
        r.NVarianceEstimates=sum(isfinite(D.VarianceMSE));
        r.VarianceBias=mean(D.VarianceBias,'omitnan');
        r.VarianceRMSE=sqrt(mean(D.VarianceMSE,'omitnan'));
        r.NormalizedVarianceRMSE=sqrt(mean(D.NormalizedVarianceMSE,'omitnan'));
        r.EstimatorVariance=mean(var([D.EstimatedVariance1 D.EstimatedVariance2 D.EstimatedVariance3],0,1,'omitnan'));
        r.EmpiricalVarianceRMSE=sqrt(mean(D.EmpiricalVarianceMSE,'omitnan'));
        r.CovarianceRelativeError=mean(D.NormalizedCovarianceError,'omitnan');
        r.NonfiniteRate=mean(D.NonfiniteEstimate); r.NegativeVarianceRate=mean(D.NegativeVariance);
        r.NFiniteCovariances=sum(~D.NonfiniteEstimate);
        if m<=4 && r.NFiniteCovariances>0
            r.IndefiniteRate=sum(D.IndefiniteCovariance)/r.NFiniteCovariances;
            [r.IndefiniteCILow,r.IndefiniteCIHigh]=wilson_ci(sum(D.IndefiniteCovariance),r.NFiniteCovariances);
        else
            r.IndefiniteRate=NaN; r.IndefiniteCILow=NaN; r.IndefiniteCIHigh=NaN;
        end
        r.EstimationFailureRate=mean(D.NonfiniteEstimate|D.NegativeVariance);
        [r.EstimationFailureCILow,r.EstimationFailureCIHigh]=wilson_ci(sum(D.NonfiniteEstimate|D.NegativeVariance),height(D));
        r.NMerge=sum(isfinite(D.MergeMSE));
        [r.MergeMSE,r.MergeMSECILow,r.MergeMSECIHigh]=mean_ci(D.MergeMSE);
        r.NormalizedMergeMSE=mean(D.NormalizedMergeMSE,'omitnan');
        r.MergeFailureRate=mean(~isfinite(D.MergeMSE));
        r.MergeR2=mean(D.MergeR2,'omitnan');
        r.MeanWeightL1=mean(D.WeightL1,'omitnan'); r.MaxAbsWeight=max(D.MaxAbsWeight);
        r.MeanUnsafeMergeMSE=mean(D.UnsafeMergeMSE,'omitnan');
        r.MaxUnsafeMergeMSE=max(D.UnsafeMergeMSE);
        r.MeanPSDThetaScale=mean(D.NCThetaScale,'omitnan');
        r.PSDRepairRate=mean(D.NCPSDRepairApplied,'omitnan');
        r.MeanNCIterations=mean(D.NCIterations,'omitnan');
        r.MeanEstimatorSeconds=mean(D.EstimatorSeconds,'omitnan');
        r.RealizedSignalErrorCorrelation=mean(D.RealizedSignalErrorCorrelation,'omitnan');
        r.TCHFittedNegativeVarianceRate=mean(D.TCHFittedNegativeVariance,'omitnan');
        r.TCHFittedMinimumEigenvalue=min(D.TCHFittedMinEigenvalue);
        r.TCHMeanReconstructionRelativeChange=mean(D.TCHReconstructionRelativeChange,'omitnan');
        if m>4
            r.NFiniteCovariances=0;
            r.NonfiniteRate=NaN; r.NegativeVarianceRate=NaN; r.IndefiniteRate=NaN;
            r.EstimationFailureRate=NaN; r.EstimationFailureCILow=NaN; r.EstimationFailureCIHigh=NaN;
            r.IndefiniteCILow=NaN; r.IndefiniteCIHigh=NaN;
        end
        if k==1
            recs=r;
        else
            recs(k)=orderfields(r,recs(1));
        end
    end
end
out=struct2table(recs);
end

function out=common_subset_comparisons(T,S)
recs=struct([]); k=0;
for s=1:height(S)
    A=T(T.Scenario==s & T.Method=="NC_PSD",:);
    for name={"TC","TCH","NC_raw"}
        B=T(T.Scenario==s & T.Method==name{1},:);
        assert(isequal(A.Replicate,B.Replicate));
        k=k+1; r=table2struct(S(s,:)); r.Contrast="NC_PSD_minus_"+name{1};
        d=A.MergeMSE-B.MergeMSE; r.NCommonMerge=sum(isfinite(d));
        [r.MergeMSEDifference,r.MergeMSEDifferenceCILow,r.MergeMSEDifferenceCIHigh]=mean_ci(d);
        d=A.VarianceMSE-B.VarianceMSE; r.NCommonVariance=sum(isfinite(d));
        [r.VarianceMSEDifference,r.VarianceMSEDifferenceCILow,r.VarianceMSEDifferenceCIHigh]=mean_ci(d);
        if k==1, recs=r; else, recs(k)=orderfields(r,recs(1)); end
    end
end
out=struct2table(recs);
end

function out=paired_temporal_summary(T,S,methods)
recs=struct([]); k=0;
for s=1:height(S)
    if S.TruthIndex(s)~=1, continue; end
    target=find(S.TruthIndex==2 & S.SNRdB==S.SNRdB(s) & S.ErrorCorrelation==S.ErrorCorrelation(s) ...
        & S.SignalErrorCorrelation==S.SignalErrorCorrelation(s));
    for m=1:numel(methods)
        A=T(T.Scenario==s & T.Method==methods{m},:);
        B=T(T.Scenario==target & T.Method==methods{m},:);
        assert(isequal(A.Replicate,B.Replicate));
        k=k+1; r=struct('SNRdB',S.SNRdB(s),'ErrorCorrelation',S.ErrorCorrelation(s), ...
            'SignalErrorCorrelation',S.SignalErrorCorrelation(s),'Method',string(methods{m}));
        d=B.MergeMSE-A.MergeMSE; r.NPairedMerge=sum(isfinite(d));
        [r.MergeMSEDifference,r.MergeMSEDifferenceCILow,r.MergeMSEDifferenceCIHigh]=mean_ci(d);
        d=B.VarianceMSE-A.VarianceMSE; r.NPairedVariance=sum(isfinite(d));
        [r.VarianceMSEDifference,r.VarianceMSEDifferenceCILow,r.VarianceMSEDifferenceCIHigh]=mean_ci(d);
        if k==1, recs=r; else, recs(k)=orderfields(r,recs(1)); end
    end
end
out=struct2table(recs);
end

function [m,lo,hi]=mean_ci(x)
x=x(isfinite(x)); m=mean(x); lo=NaN; hi=NaN;
if numel(x)>=2
    se=std(x)/sqrt(numel(x)); lo=m-1.96*se; hi=m+1.96*se;
end
end

function [lo,hi]=wilson_ci(k,n)
z=1.96; p=k/n; d=1+z^2/n; c=(p+z^2/(2*n))/d;
h=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d;
lo=max(0,c-h); hi=min(1,c+h);
end

function write_run_note(cfg,runInfo,S,P,outDir)
fid=fopen(fullfile(outDir,'RUN_NOTES.md'),'w'); assert(fid>=0);
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Executed revision synthetic experiment\n\n');
fprintf(fid,'MATLAB %s; %s; runtime %.1f s. Seed %d; %d paired Monte Carlo replicates per scenario.\n\n', ...
    runInfo.matlabVersion,runInfo.timestampUTC,runInfo.runtimeSeconds,cfg.seed,cfg.nMC);
fprintf(fid,'Three products have unit gains and heterogeneous error variances proportional to [0.5, 1, 1.5]. Each realization contains %d training and %d independently generated test days. Truth is centered and scaled to sample variance one. SNR is signal variance divided by mean target product error variance; 0.1, 3 and 10 dB equal %.5f, %.5f and %.5f in linear units. These are sensitivity settings, not empirically established typical soil-moisture SNRs.\n\n', ...
    cfg.nTrain,cfg.nTest,10^(0.1/10),10^(3/10),10^(10/10));
fprintf(fid,'Stationary truth is AR(1), phi=0.8. Seasonal truth combines a deterministic annual sinusoid (nominal variance fraction 0.5) and AR(1) processes phi=0.95 and 0.4 (fractions 0.3 and 0.2) before final normalization. This provides seasonal mean nonstationarity and two correlation scales; it is not a claim to model all soil-moisture physics or long-range dependence.\n\n');
fprintf(fid,'Base errors are stationary Gaussian AR(1), phi=0.5, with equicorrelation 0, 0.4, or 0.99 and a 1000-day burn-in. The 0.99 cases are deliberately severe near-collinearity stress tests, not representative operational cases; standard figures show 0 and 0.4, while all three are tabulated. For alpha=0, the stationary and seasonal cases use exactly the same train/test errors. For the assumption-violation scenarios, e_i=sqrt(v_i)[sqrt(1-alpha^2)u_i+alpha*y], alpha=0.3. The underlying u_i draws remain paired, but the added signal term necessarily changes realized errors. The target covariance is D[(1-alpha^2)R+alpha^2*11'']D, so total error cross-correlation is not the base correlation when alpha>0. The assumption violation is deliberately explicit.\n\n');
fprintf(fid,'All estimates use cov(X_train), never truth. Each covariance is merged using the same revised WA default PSD projection and strictly positive eigenvalue floor. NC_raw means unconstrained residual covariance BEFORE this common merge stabilization; NC_PSD additionally applies the NC decomposition-preserving rank-one shrinkage. UnsafeMergeMSE is the finite unprotected solve diagnostic, not a recommended result. Invalid estimates stay missing and are reported separately; denominators and counts appear in scenario_summary.csv.\n\n');
fprintf(fid,'TCH is exactly the existing repository TCH.m: correlated-error fminsearch objective with feasibility penalties, followed by the FIRST OUTPUT reconstruction Q-theta*theta''. It is not the simple uncorrelated pairwise least-squares hat estimator. The fourth fitted covariance output is captured for minimum-eigenvalue, negative-diagonal, and reconstruction-difference diagnostics. NC_raw uses the corrected revision optimizer before its native PSD map; it is not a bitwise reproduction of the originally submitted optimizer.\n\n');
fprintf(fid,'Performance means are conditional on successful estimation/merging; counts and failure rates must be read with them. common_subset_method_comparisons.csv additionally reports paired NC_PSD-minus-comparator differences only where both methods are valid, with the number of common replicates.\n\n');
fprintf(fid,'IndefiniteRate is conditional on finite covariance outputs (NFiniteCovariances); it is NA when every estimate is nonfinite. EstimationFailureRate and MergeFailureRate use all Monte Carlo attempts. Controls are not covariance estimators, so their covariance-admissibility rates are NA.\n\n');
fprintf(fid,'Variance metrics compare against known generating error variances; empirical covariance targets are also reported. Error-variance RMSE pools squared errors over products and Monte Carlo realizations. EstimatorVariance is the mean across products of the Monte Carlo variance of the estimated variance. Merging MSE is computed on the independent test sequence. Monte Carlo mean CIs use +/-1.96 SE across independent realizations, not across autocorrelated days; rate CIs use Wilson intervals. Paired temporal CIs use within-replicate seasonal-minus-stationary differences. These intervals quantify simulation uncertainty under this generator, not observational uncertainty.\n\n');
fprintf(fid,'## Seasonal results at 3 dB\n\n| Base corr | Signal-error corr | Method | Variance bias | Variance RMSE | Merge MSE [95%% MC CI] | Indefinite rate | Fail rate |\n|---|---|---|---|---|---|---|---|\n');
D=S(S.TruthIndex==2 & S.SNRdB==3,:);
for i=1:height(D)
    fprintf(fid,'| %.1f | %.1f | %s | %.6f | %.6f | %.6f [%.6f, %.6f] | %.1f%% | %.1f%% |\n', ...
        D.ErrorCorrelation(i),D.SignalErrorCorrelation(i),D.Method(i),D.VarianceBias(i),D.VarianceRMSE(i), ...
        D.MergeMSE(i),D.MergeMSECILow(i),D.MergeMSECIHigh(i),100*D.IndefiniteRate(i),100*D.EstimationFailureRate(i));
end
fprintf(fid,'\n## Paired seasonal minus stationary at 3 dB, uncorrelated signal/errors\n\n| Base corr | Method | Difference in merge MSE [95%% MC CI] | Difference in variance MSE [95%% MC CI] |\n|---|---|---|---|\n');
D=P(P.SNRdB==3 & P.SignalErrorCorrelation==0,:);
for i=1:height(D)
    fprintf(fid,'| %.1f | %s | %.6f [%.6f, %.6f] | %.6f [%.6f, %.6f] |\n',D.ErrorCorrelation(i),D.Method(i), ...
        D.MergeMSEDifference(i),D.MergeMSEDifferenceCILow(i),D.MergeMSEDifferenceCIHigh(i), ...
        D.VarianceMSEDifference(i),D.VarianceMSEDifferenceCILow(i),D.VarianceMSEDifferenceCIHigh(i));
end
fprintf(fid,'\nNo superiority claim should be inferred from numerical admissibility alone. Positive semidefiniteness does not establish physical identifiability, error-signal uncorrelatedness, optimality or universal ranking. The oracle and equal-weight rows are merging controls, not covariance estimators; their variance-estimation metrics are NA.\n');
end
