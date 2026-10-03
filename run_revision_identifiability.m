function out=run_revision_identifiability(nMC,outDir)
%RUN_REVISION_IDENTIFIABILITY Exact nonuniqueness and p=6 sensitivity.
% Additional controlled diagnostics; no result is an identifiability proof.
%   run_revision_identifiability(50, fullfile(pwd,'results','sensitivity'))
if nargin<1 || isempty(nMC), nMC=50; end
here=fileparts(mfilename('fullpath')); addpath(here,'-begin');
if nargin<2 || isempty(outDir), outDir=fullfile(here,'identifiability_results'); end
if ~exist(outDir,'dir'), mkdir(outDir); end
start=tic; seed=20260918;

% With p=2 these PSD decompositions have the SAME zero objective and Q.
% They have different latent loadings: unit-sum weights are shown only to
% illustrate the ambiguity in a naive downstream merge, not to assert that
% all interpretations simultaneously satisfy known unit-gain calibration.
Q2=[2 1;1 2]; c=[0.8 1 1.25]; C=cell(1,numel(c));
for k=1:numel(c)
    th=[c(k);1/c(k)]; E=Q2-th*th'; [w,wi]=WA(E);
    R=Q2-th*th'; loss=sum(abs(R(~eye(2))));
    C{k}=struct('Theta1',th(1),'Theta2',th(2),'Variance1',E(1,1), ...
        'Variance2',E(2,2),'OffDiagonalLoss',loss,'MinEigenvalue',min(eig(E)), ...
        'DecompositionResidual',norm(Q2-th*th'-E,'fro'), ...
        'UnitSumWeight1',w(1),'UnitSumWeight2',w(2),'FloorApplied',wi.regularized);
end
counterexample=struct2table([C{:}]);
writetable(counterexample,fullfile(outDir,'p2_exact_nonuniqueness.csv'));

% p=3: two PSD physical decompositions fit the same observation covariance;
% L1 selects the diagonal one although a correlated E is equally compatible
% with the data. The off-diagonal criterion is a selection assumption.
Q3=[2 1.2 1.2;1.2 2.5 1.2;1.2 1.2 3];
latent=[1 1.2]; P=cell(1,2);
for k=1:2
    th=sqrt(latent(k))*ones(3,1); E=Q3-th*th'; w=WA(E);
    P{k}=struct('SignalVariance',latent(k),'Variance1',E(1,1),'Variance2',E(2,2), ...
        'Variance3',E(3,3),'OffDiagonal12',E(1,2),'OffDiagonalLoss',sum(abs(E(~eye(3)))), ...
        'MinEigenvalue',min(eig(E)),'Weight1',w(1),'Weight2',w(2),'Weight3',w(3));
end
p3selection=struct2table([P{:}]);
writetable(p3selection,fullfile(outDir,'p3_physical_selection_ambiguity.csv'));

p=6; n=1095; burn=1000; snr=10^(3/10); rho=0.99;
v=linspace(0.5,1.5,p)'/snr;
Etrue=diag(sqrt(v))*((1-rho)*eye(p)+rho*ones(p))*diag(sqrt(v));
L=chol(Etrue,'lower');
settings=struct('Name',{'one_spectral','default_spectral','six_spectral','random8','shift_low','shift_high'}, ...
    'nStarts',{1,4,6,4,4,4},'nRandomStarts',{0,0,0,8,0,0},'betaFactor',{0.5,0.5,0.5,0.5,0.1,0.9});
floors=[1e-10 1e-8 1e-6 1e-4 1e-2]; rows=struct([]); idx=0;
for k=1:nMC
    rng(seed+k,'twister');
    y=seasonal_truth(n,burn); yt=seasonal_truth(n,burn);
    e=filter(sqrt(1-0.5^2),[1 -0.5],randn(n+burn,p)); e=e(burn+1:end,:)*L';
    et=filter(sqrt(1-0.5^2),[1 -0.5],randn(n+burn,p)); et=et(burn+1:end,:)*L';
    Q=cov(y+e); Xtest=yt+et;
    for j=1:numel(settings)
        st=settings(j);
        opts=struct('nStarts',st.nStarts,'nRandomStarts',st.nRandomStarts, ...
            'randomSeed',seed+k,'betaFactor',st.betaFactor,'enforcePSD',true);
        [E,~,~,ni]=NC(Q,opts);
        for f=floors
            [w,wi]=WA(E,struct('minEigenvalueRatio',f));
            idx=idx+1;
            row=struct('Replicate',k,'Seed',seed+k,'Setting',string(st.Name), ...
                'EigenvalueRatio',f,'nStarts',st.nStarts,'nRandomStarts',st.nRandomStarts, ...
                'BetaFactor',st.betaFactor,'Objective',ni.objective, ...
                'RawObjective',ni.objectiveBeforePSD,'PSDThetaScale',ni.thetaScalePSD, ...
                'PSDRepairApplied',ni.psdRepairApplied,'MinErrorEigenvalue',min(eig(E)), ...
                'NormalizedVarianceRMSE',sqrt(mean((diag(E)-v).^2))/mean(v), ...
                'WeightL1',sum(abs(w)),'MaxAbsWeight',max(abs(w)), ...
                'PopulationMergeMSE',w'*Etrue*w,'HeldoutMergeMSE',mean((Xtest*w-yt).^2), ...
                'FloorApplied',wi.regularized,'CovarianceRelativeChange',wi.relativeChange, ...
                'EffectiveConditionNumber',wi.effectiveConditionNumber);
            if idx==1, rows=row; else, rows(idx)=orderfields(row,rows(1)); end
        end
    end
    if mod(k,max(1,round(nMC/5)))==0 || k==1
        fprintf('Sensitivity %d/%d; %.1f seconds.\n',k,nMC,toc(start));
    end
end
replicates=struct2table(rows);
sumrows=struct([]); is=0;
for j=1:numel(settings)
    for f=floors
        D=replicates(replicates.Setting==settings(j).Name & replicates.EigenvalueRatio==f,:);
        is=is+1;
        [m,lo,hi]=mean_ci(D.HeldoutMergeMSE);
        row=struct('Setting',string(settings(j).Name),'EigenvalueRatio',f,'N',height(D), ...
            'MeanObjective',mean(D.Objective),'MeanRawObjective',mean(D.RawObjective), ...
            'MeanNormalizedVarianceRMSE',mean(D.NormalizedVarianceRMSE), ...
            'MeanWeightL1',mean(D.WeightL1),'MaxAbsWeight',max(D.MaxAbsWeight), ...
            'MeanPopulationMergeMSE',mean(D.PopulationMergeMSE), ...
            'MeanHeldoutMergeMSE',m,'HeldoutMSECILow',lo,'HeldoutMSECIHigh',hi, ...
            'PSDRepairRate',mean(D.PSDRepairApplied),'FloorRate',mean(D.FloorApplied), ...
            'MeanCovarianceRelativeChange',mean(D.CovarianceRelativeChange));
        if is==1, sumrows=row; else, sumrows(is)=orderfields(row,sumrows(1)); end
    end
end
summary=struct2table(sumrows);
cfg=struct('seed',seed,'nMC',nMC,'p',p,'nTrain',n,'nTest',n,'burn',burn,'snrDb',3,'rho',rho, ...
    'floors',floors,'settings',settings,'runtimeSeconds',toc(start),'matlabVersion',version);
writetable(replicates,fullfile(outDir,'multistart_floor_replicates.csv'));
writetable(summary,fullfile(outDir,'multistart_floor_summary.csv'));
save(fullfile(outDir,'identifiability_sensitivity.mat'),'counterexample','p3selection','replicates','summary','cfg','Q2','Q3','Etrue');
fid=fopen(fullfile(outDir,'RUN_NOTES.md'),'w'); assert(fid>=0); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Identifiability and numerical sensitivity diagnostics\n\n');
fprintf(fid,'MATLAB %s; seed %d; %d independent realizations; runtime %.1f s.\n\n',version,seed,nMC,cfg.runtimeSeconds);
fprintf(fid,'The exact p=2 example uses Q=[2,1;1,2] and theta=[c,1/c], c=0.8,1,1.25. All three residuals are positive definite and have zero off-diagonal objective, but different residual variances and unit-sum weights. This establishes nonuniqueness even after PSD enforcement. The candidate theta vectors correspond to different loadings; the displayed unit-sum weights illustrate downstream ambiguity when scaling assumptions are not fixed. Known equal loadings would impose additional information and eliminate this particular family.\n\n');
fprintf(fid,'The p=3 example uses Q=[2,1.2,1.2;1.2,2.5,1.2;1.2,1.2,3]. Both signal variance 1 with correlated residual errors and signal variance 1.2 with diagonal errors generate exactly Q and are PSD. The off-diagonal loss favors the diagonal error interpretation. Therefore, numerical selection by the NC loss does not establish unique physical recovery, even for p>=3.\n\n');
fprintf(fid,'The p=6 finite-sample sensitivity uses unit gains, seasonal truth (annual sine + AR 0.95/0.4), signal-uncorrelated AR(1) errors phi=0.5, error cross-correlation 0.99, heterogeneous variances proportional to [0.5,...,1.5], and SNR 3 dB. The near-perfect correlation is a deliberately severe stress setting intended to expose PSD-boundary and eigenvalue-floor effects, not an operationally representative claim. Training and independent test lengths are each 1095. Six initialization settings and five relative eigenvalue floors are applied to the same covariance per realization. Variations quantify algorithmic sensitivity under this controlled setup; stable numerical results do not prove physical identifiability.\n\n');
fprintf(fid,'| Initialization | Floor ratio | Mean loss | Heldout MSE [95%% MC CI] | PSD repair rate | Mean weight L1 |\n|---|---|---|---|---|---|\n');
for j=1:height(summary)
    fprintf(fid,'| %s | %.0e | %.6g | %.6g [%.6g, %.6g] | %.1f%% | %.6g |\n',summary.Setting(j),summary.EigenvalueRatio(j), ...
        summary.MeanObjective(j),summary.MeanHeldoutMergeMSE(j),summary.HeldoutMSECILow(j),summary.HeldoutMSECIHigh(j), ...
        100*summary.PSDRepairRate(j),summary.MeanWeightL1(j));
end
out=struct('summary',summary,'counterexample',counterexample,'cfg',cfg,'outDir',outDir);
fprintf('Sensitivity complete: %s\n',outDir);
end

function y=seasonal_truth(n,burn)
z=randn(n+burn,2);
s=filter(sqrt(1-0.95^2),[1 -0.95],z(:,1)); s=s(burn+1:end);
f=filter(sqrt(1-0.4^2),[1 -0.4],z(:,2)); f=f(burn+1:end);
y=sqrt(0.5)*sqrt(2)*sin(2*pi*(1:n)'/365)+sqrt(0.3)*s+sqrt(0.2)*f;
y=(y-mean(y))/std(y);
end

function [m,lo,hi]=mean_ci(x)
x=x(isfinite(x)); m=mean(x); se=std(x)/sqrt(numel(x)); lo=m-1.96*se; hi=m+1.96*se;
end
