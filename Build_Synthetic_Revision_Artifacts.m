function report = Build_Synthetic_Revision_Artifacts(mainDir, sensitivityDir, artifactDir)
%BUILD_SYNTHETIC_REVISION_ARTIFACTS Tables, publication figures and text.
% report = Build_Synthetic_Revision_Artifacts(mainDir,sensitivityDir,artifactDir)
% Reads completed revision CSVs only; never reruns or changes an estimator.
% artifactDir must not exist or must be empty. Requires base MATLAB.
% PNG figures are 600 dpi; matching PDFs use vector output. Every statistic,
% confidence interval and sample count is read from the supplied CSVs.

arguments
    mainDir {mustBeTextScalar}
    sensitivityDir {mustBeTextScalar}
    artifactDir {mustBeTextScalar}
end
mainDir=char(mainDir); sensitivityDir=char(sensitivityDir); artifactDir=char(artifactDir);
assert(isfolder(mainDir),'BuildSynthetic:MissingMain','Main results directory does not exist.');
assert(isfolder(sensitivityDir),'BuildSynthetic:MissingSensitivity','Sensitivity results directory does not exist.');
if isfolder(artifactDir)
    d=dir(artifactDir); d=d(~ismember({d.name},{'.','..'}));
    assert(isempty(d),'BuildSynthetic:OutputNotEmpty', ...
        'artifactDir must be empty; existing files will not be overwritten: %s',artifactDir);
else
    [ok,msg]=mkdir(artifactDir); assert(ok,'BuildSynthetic:CreateOutput','%s',msg);
end

S=read_required(mainDir,'scenario_summary.csv');
P=read_required(mainDir,'paired_seasonal_minus_stationary.csv');
C=read_required(mainDir,'common_subset_method_comparisons.csv');
D=read_required(mainDir,'scenario_configuration.csv');
I=read_required(sensitivityDir,'multistart_floor_summary.csv');
P2=read_required(sensitivityDir,'p2_exact_nonuniqueness.csv');
P3=read_required(sensitivityDir,'p3_physical_selection_ambiguity.csv');
cfg=struct(); meta=struct();
cfgFile=fullfile(mainDir,'run_configuration.json');
if isfile(cfgFile)
    J=jsondecode(fileread(cfgFile));
    if isfield(J,'config'), cfg=J.config; end
    if isfield(J,'runInfo'), meta=J.runInfo; end
end
validate_columns(S,{'Scenario','Truth','SNRdB','ErrorCorrelation','SignalErrorCorrelation','Method', ...
    'NReplicates','NVarianceEstimates','NFiniteCovariances','NMerge','VarianceBias','VarianceRMSE', ...
    'NormalizedVarianceRMSE','EstimationFailureRate','IndefiniteRate','MergeFailureRate', ...
    'MergeMSE','MergeMSECILow','MergeMSECIHigh','PSDRepairRate','MeanPSDThetaScale'},'scenario_summary');
validate_columns(P,{'SNRdB','ErrorCorrelation','SignalErrorCorrelation','Method','NPairedMerge', ...
    'MergeMSEDifference','MergeMSEDifferenceCILow','MergeMSEDifferenceCIHigh'},'paired temporal');
validate_columns(C,{'Truth','SNRdB','ErrorCorrelation','SignalErrorCorrelation','Contrast','NCommonMerge', ...
    'MergeMSEDifference','MergeMSEDifferenceCILow','MergeMSEDifferenceCIHigh'},'common subset');
validate_columns(I,{'Setting','EigenvalueRatio','N','MeanObjective','MeanRawObjective','MeanHeldoutMergeMSE', ...
    'HeldoutMSECILow','HeldoutMSECIHigh','MeanWeightL1','PSDRepairRate','FloorRate'},'sensitivity');
S.Method=string(S.Method); S.Truth=string(S.Truth); P.Method=string(P.Method);
C.Contrast=string(C.Contrast); C.Truth=string(C.Truth); I.Setting=string(I.Setting);

methods=["TC","TCH","NC_raw","NC_PSD"];
stress=S(S.ErrorCorrelation>=0.9 & ismember(S.Method,methods),:);
files={};
tableSpecs={D,'Table_R1_Experimental_Design.csv'; S,'Table_R2_Synthetic_Performance.csv'; ...
    P,'Table_R3_Paired_Temporal_Changes.csv'; C,'Table_R4_Common_Subset_Contrasts.csv'; ...
    stress,'Table_R5_PSD_Stress.csv'; I,'Table_R6_Initialization_Floor_Sensitivity.csv'; ...
    P2,'Table_R7_P2_Nonuniqueness.csv'; P3,'Table_R8_P3_Physical_Ambiguity.csv'};
for k=1:size(tableSpecs,1)
    path=fullfile(artifactDir,tableSpecs{k,2}); writetable(tableSpecs{k,1},path); files{end+1}=path; %#ok<AGROW>
end

alphaValues=finite_unique(S.SignalErrorCorrelation);
if isempty(alphaValues), alphaValues=NaN; end
normalRho=finite_unique(S.ErrorCorrelation(S.ErrorCorrelation<0.9));
if isempty(normalRho), normalRho=finite_unique(S.ErrorCorrelation); end
truths=unique(S.Truth,'stable');
if isempty(truths), truths="NA"; end
if isempty(normalRho), normalRho=NaN; end
for a=alphaValues(:)'
    f=variance_figure(S,truths,normalRho,a,["TC","TCH","NC_PSD"]);
    files=[files,save_figure(f,artifactDir,['Figure_R1_Variance_RMSE_alpha_' numeric_tag(a)])]; %#ok<AGROW>
end
truth=preferred_truth(truths);
f=merge_figure(S,truth,normalRho,alphaValues,["TC","TCH","NC_PSD"]);
files=[files,save_figure(f,artifactDir,'Figure_R2_Heldout_Merging')];
f=stress_figure(S,C,truth,methods);
files=[files,save_figure(f,artifactDir,'Figure_R3_Stress_Feasibility_and_PSD_Cost')];
[f,initSetting,referenceFloor]=sensitivity_figure(I);
files=[files,save_figure(f,artifactDir,'Figure_R4_Initialization_and_Eigenfloor')];
f=initialization_figure(I,referenceFloor);
files=[files,save_figure(f,artifactDir,'Figure_R4b_Initialization')];
f=merge_figure(S,truth,normalRho,alphaValues,["NC_PSD","equal_weight","oracle"]);
files=[files,save_figure(f,artifactDir,'Figure_R5_Merging_Controls')];
f=temporal_figure(P,normalRho,alphaValues);
files=[files,save_figure(f,artifactDir,'Figure_R6_Paired_Seasonal_Change')];
files=[files,Draw_JoH_Graphical_Abstract(artifactDir)];

counts=struct('mainScenarioCount',numel(unique(S.Scenario)), ...
    'mainEstimatorAttempts',sum(S.NReplicates(ismember(S.Method,methods))), ...
    'mainSummaryRows',height(S),'pairedTemporalRows',height(P), ...
    'commonSubsetRows',height(C),'sensitivitySummaryRows',height(I));
reportPath=fullfile(artifactDir,'synthetic_results.md');
write_results(reportPath,S,P,C,I,P2,P3,cfg,meta,counts,truth,initSetting,referenceFloor);
files{end+1}=reportPath;
guidePath=fullfile(artifactDir,'ARTIFACT_README.md');
write_guide(guidePath,files,cfg);
files{end+1}=guidePath;
report=struct('files',{files},'counts',counts,'mainDir',mainDir, ...
    'sensitivityDir',sensitivityDir,'artifactDir',artifactDir, ...
    'initializationReferenceSetting',char(initSetting),'initializationReferenceFloor',referenceFloor);
manifestPath=fullfile(artifactDir,'artifact_manifest.json');
report.files{end+1}=manifestPath;
fid=fopen(manifestPath,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('Synthetic artifacts: %d files written to %s\n',numel(report.files),artifactDir);
end

function T=read_required(folder,name)
path=fullfile(folder,name);
assert(isfile(path),'BuildSynthetic:MissingCSV','Required completed-run CSV is missing: %s',path);
T=readtable(path,'TextType','string','VariableNamingRule','preserve');
end

function validate_columns(T,names,label)
missing=setdiff(names,T.Properties.VariableNames);
assert(isempty(missing),'BuildSynthetic:MissingColumns','%s missing columns: %s',label,strjoin(missing,', '));
end

function u=finite_unique(v)
u=unique(v(isfinite(v))); u=sort(u(:));
end

function tag=numeric_tag(v)
if isnan(v), tag='NA'; else, tag=strrep(strrep(sprintf('%.6g',v),'-','minus'),'.','p'); end
end

function truth=preferred_truth(truths)
idx=find(contains(lower(truths),'seasonal'),1);
if isempty(idx), idx=1; end
truth=truths(idx);
end

function f=new_figure(w,h)
f=JoHFigureStyle('figure',w,h);
end

function f=variance_figure(S,truths,rhos,alpha,methods)
nr=numel(truths); nc=numel(rhos);
f=new_figure(max(20,10*nc),max(11,8*nr));
tl=tiledlayout(f,nr,nc,'Padding','compact','TileSpacing','compact');
for t=1:nr
    for r=1:nc
        ax=nexttile(tl); hold(ax,'on'); present=false;
        for m=1:numel(methods)
            B=S(S.Truth==truths(t) & eqnum(S.ErrorCorrelation,rhos(r)) & ...
                eqnum(S.SignalErrorCorrelation,alpha) & S.Method==methods(m),:);
            B=sortrows(B,'SNRdB');
            if isempty(B) || ~any(isfinite(B.NormalizedVarianceRMSE)), continue; end
            plot(ax,B.SNRdB,B.NormalizedVarianceRMSE,style(methods(m)), ...
                'Color',method_color(methods(m)),'LineWidth',1.4,'MarkerSize',6, ...
                'DisplayName',method_label(methods(m),false)); present=true;
        end
        if ~present, show_na(ax,'NA: no matching covariance estimates'); end
        title(ax,sprintf('(%c) %s | rho = %s',96+(t-1)*nc+r,display_truth(truths(t)),num(rhos(r))));
        xlabel(ax,'SNR (dB)'); ylabel(ax,'Normalized variance RMSE'); axes_style(ax);
        if present, pad_x(ax,S.SNRdB); end
        if t==1 && r==1 && present, JoHFigureStyle('legend',ax,tl); end
    end
end
title(tl,sprintf('Signal-error parameter alpha = %s',num(alpha)),'FontSize',13,'FontWeight','normal');
end

function f=merge_figure(S,truth,rhos,alphas,methods)
f=new_figure(max(20,10*numel(rhos)),max(11,8*numel(alphas)));
tl=tiledlayout(f,numel(alphas),numel(rhos),'Padding','compact','TileSpacing','compact');
allMethods=methods;
for a=1:numel(alphas)
    for r=1:numel(rhos)
        ax=nexttile(tl); hold(ax,'on'); present=false;
        for m=1:numel(allMethods)
            B=S(S.Truth==truth & eqnum(S.ErrorCorrelation,rhos(r)) & ...
                eqnum(S.SignalErrorCorrelation,alphas(a)) & S.Method==allMethods(m),:);
            B=sortrows(B,'SNRdB'); if isempty(B) || ~any(isfinite(B.MergeMSE)), continue; end
            draw_errorbar(ax,B.SNRdB,B.MergeMSE,B.MergeMSECILow,B.MergeMSECIHigh, ...
                method_color(allMethods(m)),style(allMethods(m)),method_label(allMethods(m),true));
            present=true;
        end
        if ~present, show_na(ax,'NA: no matching test-set estimates'); end
        title(ax,sprintf('(%c) rho = %s | alpha = %s',96+(a-1)*numel(rhos)+r,num(rhos(r)),num(alphas(a))));
        xlabel(ax,'SNR (dB)'); ylabel(ax,'Test merging MSE'); axes_style(ax);
        if present, pad_x(ax,S.SNRdB); end
        if a==1 && r==1 && present, JoHFigureStyle('legend',ax,tl); end
    end
end
title(tl,[display_truth(truth) ' truth'],'FontSize',13,'FontWeight','normal','Interpreter','none');
end

function f=stress_figure(S,C,truth,methods)
f=new_figure(22,18); tl=tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
rhos=finite_unique(S.ErrorCorrelation(S.ErrorCorrelation>=0.9));
if isempty(rhos)
    for k=1:4, ax=nexttile(tl); show_na(ax,'NA: no high-correlation stress scenarios'); end
    title(tl,'Stress diagnostics unavailable'); return;
end
rho=max(rhos); alpha=minfinite(S.SignalErrorCorrelation,NaN);
B=S(S.Truth==truth & eqnum(S.ErrorCorrelation,rho) & eqnum(S.SignalErrorCorrelation,alpha),:);
rateCols={'EstimationFailureRate','IndefiniteRate'};
for panel=1:2
    ax=nexttile(tl); hold(ax,'on'); missingLabels=strings(0,1);
    for m=1:numel(methods)
        A=sortrows(B(B.Method==methods(m),:),'SNRdB');
        if isempty(A) || ~any(isfinite(A.(rateCols{panel})))
            missingLabels(end+1)=method_label(methods(m),false); %#ok<AGROW>
            continue;
        end
        plot(ax,A.SNRdB,100*A.(rateCols{panel}),style(methods(m)), ...
            'Color',method_color(methods(m)),'LineWidth',1.5,'MarkerSize',6, ...
            'DisplayName',method_label(methods(m),false));
    end
    ylim(ax,[-3 108]); xlabel(ax,'SNR (dB)'); ylabel(ax,'Rate (%)'); axes_style(ax);
    pad_x(ax,B.SNRdB);
    if panel==1, title(ax,'(a) Failed estimation / attempts');
    else, title(ax,'(b) Indefinite / finite outputs'); end
    if panel==1, JoHFigureStyle('legend',ax,tl); end
    if ~isempty(missingLabels)
        text(ax,0.02,0.04,"NA: "+strjoin(missingLabels,', '),'Units','normalized', ...
            'Interpreter','none','FontSize',11,'BackgroundColor','w');
    end
end
ax=nexttile(tl); hold(ax,'on');
A=sortrows(C(C.Truth==truth & eqnum(C.ErrorCorrelation,rho) & ...
    eqnum(C.SignalErrorCorrelation,alpha) & C.Contrast=="NC_PSD_minus_NC_raw",:),'SNRdB');
if isempty(A) || ~any(A.NCommonMerge>0)
    show_na(ax,'NA: no common valid raw/PSD merges');
else
    draw_errorbar(ax,A.SNRdB,A.MergeMSEDifference,A.MergeMSEDifferenceCILow,A.MergeMSEDifferenceCIHigh, ...
        method_color('NC_PSD'),'-d','NC PSD - NC raw'); yline(ax,0,':','HandleVisibility','off');
end
title(ax,'(c) PSD map: paired MSE change');
xlabel(ax,'SNR (dB)'); ylabel(ax,'MSE difference (PSD - raw)'); axes_style(ax);
if ~isempty(A), pad_x(ax,A.SNRdB); end
ax=nexttile(tl); hold(ax,'on');
A=sortrows(B(B.Method=="NC_PSD",:),'SNRdB');
if isempty(A)
    show_na(ax,'NA: no native-PSD diagnostics');
else
    plot(ax,A.SNRdB,A.MeanPSDThetaScale,'-s','LineWidth',1.5,'Color',method_color('NC_PSD'));
    ylabel(ax,'Mean signal-vector scale');
    pad_x(ax,A.SNRdB);
end
title(ax,'(d) Size of native PSD correction'); xlabel(ax,'SNR (dB)'); axes_style(ax);
title(tl,sprintf('Stress: rho = %s | alpha = %s',num(rho),num(alpha)),'FontSize',13,'FontWeight','normal');
end

function [f,setting,refFloor]=sensitivity_figure(I)
f=new_figure(20,10); tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
settings=unique(I.Setting,'stable'); floors=finite_unique(I.EigenvalueRatio);
assert(~isempty(settings)&&~isempty(floors),'No sensitivity results.');
if any(settings=="default_spectral"), setting="default_spectral"; else, setting=settings(1); end
[~,j]=min(abs(log10(floors)-log10(1e-8))); refFloor=floors(j);
A=sortrows(I(I.Setting==setting,:),'EigenvalueRatio');
ax=nexttile(tl); hold(ax,'on');
draw_errorbar(ax,A.EigenvalueRatio,A.MeanHeldoutMergeMSE,A.HeldoutMSECILow,A.HeldoutMSECIHigh, ...
    method_color('NC_PSD'),'-o','NC'); set(ax,'XScale','log'); pad_x(ax,A.EigenvalueRatio);
title(ax,'(a) Prediction'); xlabel(ax,'Relative eigenvalue floor'); ylabel(ax,'Test merging MSE'); axes_style(ax);
ax=nexttile(tl); plot(ax,A.EigenvalueRatio,A.MeanWeightL1,'-o','LineWidth',1.8,'Color',method_color('NC_PSD'));
set(ax,'XScale','log'); title(ax,'(b) Weight magnitude'); xlabel(ax,'Relative eigenvalue floor');
ylabel(ax,'Mean sum of absolute weights'); axes_style(ax); pad_x(ax,A.EigenvalueRatio);
end

function f=initialization_figure(I,refFloor)
f=new_figure(22,12); tl=tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
B=I(eqnum(I.EigenvalueRatio,refFloor),:); n=height(B); x=(1:n)';
labels=strrep(strrep(strrep(B.Setting,'default_spectral','Spectral'),'_',' '),'perturbation','perturb.');
ax=nexttile(tl); hold(ax,'on');
for j=1:n
    plot(ax,[B.HeldoutMSECILow(j) B.HeldoutMSECIHigh(j)],[j j],'-','Color',method_color('NC_PSD'),'LineWidth',1.8);
end
plot(ax,B.MeanHeldoutMergeMSE,x,'o','Color',method_color('NC_PSD'),'MarkerFaceColor',method_color('NC_PSD'));
set(ax,'YTick',x,'YTickLabel',labels,'TickLabelInterpreter','none','YDir','reverse'); ylim(ax,[.4 n+.6]);
title(ax,'(a) Prediction'); xlabel(ax,'Test merging MSE'); axes_style(ax);
ax=nexttile(tl); hold(ax,'on');
plot(ax,B.MeanRawObjective,x-.10,'o','Color',method_color('NC_raw'),'MarkerFaceColor',method_color('NC_raw'),'DisplayName','Before PSD');
plot(ax,B.MeanObjective,x+.10,'d','Color',method_color('NC_PSD'),'MarkerFaceColor',method_color('NC_PSD'),'DisplayName','After PSD');
set(ax,'YTick',x,'YTickLabel',[],'YDir','reverse'); ylim(ax,[.4 n+.6]);
title(ax,'(b) Covariance fit'); xlabel(ax,'Mean off-diagonal loss'); axes_style(ax);
JoHFigureStyle('legend',ax,tl);
end

function f=temporal_figure(P,rhos,alphas)
f=new_figure(20,16); tl=tiledlayout(f,numel(alphas),numel(rhos),'Padding','compact','TileSpacing','compact');
methods=["TC","TCH","NC_PSD"];
for a=1:numel(alphas)
    for r=1:numel(rhos)
        ax=nexttile(tl); hold(ax,'on');
        for m=1:numel(methods)
            B=sortrows(P(eqnum(P.ErrorCorrelation,rhos(r)) & eqnum(P.SignalErrorCorrelation,alphas(a)) & P.Method==methods(m),:),'SNRdB');
            draw_errorbar(ax,B.SNRdB,B.MergeMSEDifference,B.MergeMSEDifferenceCILow,B.MergeMSEDifferenceCIHigh, ...
                method_color(methods(m)),style(methods(m)),method_label(methods(m),false));
        end
        yline(ax,0,':','HandleVisibility','off'); axes_style(ax); pad_x(ax,P.SNRdB);
        xlabel(ax,'SNR (dB)'); ylabel(ax,'Merging MSE change');
        title(ax,sprintf('(%c) rho = %s | alpha = %s',96+(a-1)*numel(rhos)+r,num(rhos(r)),num(alphas(a))));
        if a==1&&r==1, JoHFigureStyle('legend',ax,tl); end
    end
end
end

function draw_errorbar(ax,x,y,lo,hi,color,lineSpec,label)
ok=isfinite(x) & isfinite(y);
if ~any(ok), return; end
x=x(ok); y=y(ok); lo=lo(ok); hi=hi(ok);
down=max(0,y-lo); up=max(0,hi-y);
% An unavailable CI is displayed as a point with no error bar, not as zero uncertainty in text.
down(~isfinite(down))=0; up(~isfinite(up))=0;
errorbar(ax,x,y,down,up,lineSpec,'Color',color,'LineWidth',1.8, ...
    'MarkerSize',6,'CapSize',4,'DisplayName',label);
end

function axes_style(ax)
JoHFigureStyle('axes',ax);
end

function pad_x(ax,x)
x=x(isfinite(x)); if isempty(x), return; end
if strcmp(ax.XScale,'log')
    x=x(x>0); if isempty(x), return; end
    low=log10(min(x)); high=log10(max(x)); pad=max(.04*(high-low),.08);
    xlim(ax,10.^[low-pad,high+pad]);
else
    low=min(x); high=max(x); pad=max(.03*(high-low),.05*max(abs([low high])));
    if pad==0, pad=.1; end
    xlim(ax,[low-pad,high+pad]);
end
end

function show_na(ax,label)
axis(ax,[0 1 0 1]); axis(ax,'off'); text(ax,.5,.5,label,'HorizontalAlignment','center', ...
    'Interpreter','none','FontSize',9,'Color',[.35 .35 .35]);
end

function c=method_color(name)
c=JoHFigureStyle('color',name);
end

function v=style(name)
switch string(name)
    case "NC_raw", v='--^';
    case "NC_PSD", v='-d';
    case "equal_weight", v='--s';
    case "oracle", v=':x';
    case "TC", v='-o';
    case "TCH", v='--s';
    otherwise, v='-o';
end
end

function s=method_label(name,merging)
name=string(name);
switch name
    case "NC_raw", s='NC raw';
    case "NC_PSD", s='NC PSD';
    case "TCH", s='TCH';
    case "equal_weight", s='Equal weights';
    case "oracle", s='Oracle covariance';
    otherwise, s=char(name);
end
% The common SPD safeguard and TCH variant are documented in the caption.
end

function s=display_truth(v)
if contains(string(v),'seasonal'), s='Seasonal'; elseif contains(string(v),'stationary'), s='Stationary'; else, s=char(strrep(string(v),'_',' ')); end
end

function tf=eqnum(x,v)
if ~isfinite(v), tf=false(size(x)); else, tf=abs(x-v)<=32*eps(max(1,abs(v))); end
end

function paths=save_figure(f,folder,name)
clean=onCleanup(@()close(f)); %#ok<NASGU>
paths=JoHFigureStyle('export',f,folder,name);
end

function v=minfinite(x,fallback)
x=x(isfinite(x)); if isempty(x), v=fallback; else, v=min(x); end
end

function s=num(x)
if isempty(x) || ~isscalar(x) || ~isfinite(x), s='NA'; else, s=sprintf('%.6g',x); end
end

function s=ci(m,lo,hi)
s=sprintf('%s [%s, %s]',num(m),num(lo),num(hi));
end

function s=count_range(v)
v=v(isfinite(v)); if isempty(v), s='NA'; elseif min(v)==max(v), s=num(v(1));
else, s=[num(min(v)) ' to ' num(max(v))]; end
end

function write_results(path,S,P,C,I,P2,P3,cfg,meta,counts,truth,setting,refFloor)
fid=fopen(path,'w'); assert(fid>=0); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Executed synthetic experiment results\n\n');
fprintf(fid,'This report was generated from completed CSV outputs. It contains no global or ISMN results. The experiment contains **%s scenarios**, **%s Monte Carlo attempts per scenario/method**, and **%s covariance-estimator attempts** in total. The sensitivity summary contains %s rows with %s realizations per setting.\n\n', ...
    num(counts.mainScenarioCount),count_range(S.NReplicates),num(counts.mainEstimatorAttempts), ...
    num(height(I)),count_range(I.N));
fprintf(fid,'All main merging results apply the same WA SPD safeguard to TC, the repository TCH variant, NC raw, and NC PSD. NC raw denotes the current optimizer before the native PSD feasibility map; it is not an unprotected matrix inversion. TCH denotes the repository penalty optimizer''s first reconstructed covariance output, not its fourth fitted covariance output.\n\n');
if isfield(cfg,'nTrain') && isfield(cfg,'nTest')
    fprintf(fid,'Each realization uses %s training samples and %s independently generated test samples. The recorded seed is %s. ',num(cfg.nTrain),num(cfg.nTest),fieldnum(cfg,'seed'));
end
if isfield(cfg,'SNRdB')
    fprintf(fid,'The configured SNR values are %s dB, corresponding to signal/mean-error-power ratios %s. ', ...
        strjoin(compose('%.6g',cfg.SNRdB(:)'),', '),strjoin(compose('%.6g',10.^(cfg.SNRdB(:)'/10)),', '));
end
fprintf(fid,'These SNR settings are sensitivity scenarios, not claims about typical operational retrieval SNR.\n\n');
if isfield(meta,'exportRecoveredFromCSV') && meta.exportRecoveredFromCSV
    fprintf(fid,'The source archive records recovery of artifact export from already completed CSV metrics; numerical results retain CSV precision.\n\n');
end
fprintf(fid,'## Denominators and interpretation\n\n');
fprintf(fid,'Estimation and merge failure rates use every attempted realization. Indefiniteness rates use only finite covariance outputs (`NFiniteCovariances`), and are NA if that denominator is zero. Accuracy means are conditional on available estimates and are accompanied by valid counts. Equal-weight and known-covariance controls are not covariance estimators; covariance-estimation metrics for them are NA. The paired comparison table uses only realizations where both compared methods are valid. Confidence intervals quantify Monte Carlo uncertainty across independent realizations, not independent daily samples.\n\n');
fprintf(fid,'## Paired temporal comparison\n\n');
fprintf(fid,'Positive differences below mean larger seasonal than stationary error. Read the CI and common-valid count together; a small change does not establish universal immunity to seasonality.\n\n');
fprintf(fid,'| SNR dB | Base error corr. | Signal-error corr. | Method | Seasonal minus stationary test MSE [95%% CI] | Common merge n |\n|---|---|---|---|---|---|\n');
if isempty(P), fprintf(fid,'| NA | NA | NA | NA | NA | 0 |\n'); end
for k=1:height(P)
    if ~ismember(P.Method(k),["TC","TCH","NC_PSD"]), continue; end
    fprintf(fid,'| %s | %s | %s | %s | %s | %s |\n',num(P.SNRdB(k)),num(P.ErrorCorrelation(k)), ...
        num(P.SignalErrorCorrelation(k)),char(P.Method(k)),ci(P.MergeMSEDifference(k),P.MergeMSEDifferenceCILow(k),P.MergeMSEDifferenceCIHigh(k)),num(P.NPairedMerge(k)));
end
fprintf(fid,'\n## Primary NC versus comparator on common valid realizations\n\n');
fprintf(fid,'These are NC PSD minus comparator differences. Negative values favor NC PSD for that metric and subset; positive values favor the comparator. They do not compare successful NC estimates with failed comparator attempts.\n\n');
snrs=finite_unique(C.SNRdB);
if ~isempty(snrs)
    representativeSNR=snrs(ceil(numel(snrs)/2));
    representativeAlpha=minfinite(C.SignalErrorCorrelation,NaN);
    R=C(C.Truth==truth & eqnum(C.SNRdB,representativeSNR) & ...
        eqnum(C.SignalErrorCorrelation,representativeAlpha) & ...
        ismember(C.Contrast,["NC_PSD_minus_TC","NC_PSD_minus_TCH"]),:);
    fprintf(fid,'At the middle tested SNR (%s dB), using %s truth and signal-error correlation %s:\n\n', ...
        num(representativeSNR),display_truth(truth),num(representativeAlpha));
    for k=1:height(R)
        r=R(k,:);
        fprintf(fid,'- At base error correlation %s, `%s` gives a test-MSE difference of %s on %s common valid realizations.\n', ...
            num(r.ErrorCorrelation),char(r.Contrast),ci(r.MergeMSEDifference,r.MergeMSEDifferenceCILow,r.MergeMSEDifferenceCIHigh),num(r.NCommonMerge));
    end
    if isempty(R), fprintf(fid,'- NA: there are no matching common-subset contrasts.\n'); end
    fprintf(fid,'\n');
end
fprintf(fid,'| Truth | SNR dB | Base corr. | Signal-error corr. | Contrast | Test MSE difference [95%% CI] | Common n |\n|---|---|---|---|---|---|---|\n');
selected=C(C.Truth==truth,:);
if isempty(selected), fprintf(fid,'| NA | NA | NA | NA | NA | NA | 0 |\n'); end
for k=1:height(selected)
    r=selected(k,:);
    fprintf(fid,'| %s | %s | %s | %s | %s | %s | %s |\n',char(r.Truth),num(r.SNRdB),num(r.ErrorCorrelation), ...
        num(r.SignalErrorCorrelation),char(r.Contrast),ci(r.MergeMSEDifference,r.MergeMSEDifferenceCILow,r.MergeMSEDifferenceCIHigh),num(r.NCommonMerge));
end
fprintf(fid,'\n## Severe error-correlation stress and PSD cost\n\n');
fprintf(fid,'Error correlations at or above 0.9 are deliberately severe stress settings. The native PSD map and the WA SPD floor are distinct steps. In the raw-versus-PSD contrast both estimates pass through the same WA safeguard; its MSE difference therefore measures the additional native PSD map, not a comparison against an unprotected inverse. Numerical admissibility does not imply recovery of the physical generating covariance or improved merging accuracy.\n\n');
fprintf(fid,'| Truth | SNR | Base corr. | Signal-error corr. | Method | Attempts | Finite covariances | Estimation failure %% | Indefinite %% | PSD repair %% | Mean theta scale | Test MSE [95%% CI] | Valid merge n |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
T=S(S.Truth==truth & S.ErrorCorrelation>=0.9 & ismember(S.Method,["TC","TCH","NC_raw","NC_PSD"]),:);
if isempty(T), fprintf(fid,'| NA | NA | NA | NA | NA | 0 | 0 | NA | NA | NA | NA | NA | 0 |\n'); end
for k=1:height(T)
    r=T(k,:);
    fprintf(fid,'| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
        char(r.Truth),num(r.SNRdB),num(r.ErrorCorrelation),num(r.SignalErrorCorrelation),char(r.Method),num(r.NReplicates), ...
        num(r.NFiniteCovariances),num(100*r.EstimationFailureRate),num(100*r.IndefiniteRate),num(100*r.PSDRepairRate), ...
        num(r.MeanPSDThetaScale),ci(r.MergeMSE,r.MergeMSECILow,r.MergeMSECIHigh),num(r.NMerge));
end
fprintf(fid,'\n## Initialization and eigenvalue-floor sensitivity\n\n');
fprintf(fid,'The figure uses initialization `%s` for its floor sweep and floor `%s` for its initialization comparison. The table retains every tested combination. Stability of this numerical experiment does not resolve physical nonidentifiability. The two-stage PSD map can increase the off-diagonal loss, so more starts do not guarantee a better post-map objective or better prediction.\n\n',char(setting),num(refFloor));
fprintf(fid,'| Setting | Relative floor | n | Mean loss before/after PSD | Test MSE [95%% CI] | PSD repair %% | Floor %% | Mean weight L1 |\n|---|---|---|---|---|---|---|---|\n');
if isempty(I), fprintf(fid,'| NA | NA | 0 | NA | NA | NA | NA | NA |\n'); end
for k=1:height(I)
    r=I(k,:); fprintf(fid,'| %s | %s | %s | %s / %s | %s | %s | %s | %s |\n',char(r.Setting),num(r.EigenvalueRatio),num(r.N), ...
        num(r.MeanRawObjective),num(r.MeanObjective),ci(r.MeanHeldoutMergeMSE,r.HeldoutMSECILow,r.HeldoutMSECIHigh), ...
        num(100*r.PSDRepairRate),num(100*r.FloorRate),num(r.MeanWeightL1));
end
fprintf(fid,'\n## Identifiability and numerical-optimization qualifications\n\n');
fprintf(fid,'The exported p=2 and p=3 tables contain %s and %s candidate decompositions, respectively. The p=2 example permits distinct PSD residuals with the same loss; its candidates imply different loadings, so known equal gains would supply additional identifying information. The p=3 common-mode example demonstrates covariance ambiguity while its unit-sum optimal weights remain unchanged. Thus nonidentifiability does not force weight changes in every special case.\n\n',num(height(P2)),num(height(P3)));
fprintf(fid,'For a compatible p=3 covariance with nonnegative TC residual variances, the TC signal products can supply a zero-off-diagonal-loss feasible candidate. A positive final NC loss may therefore include numerical optimization approximation. The projected-subgradient stopping rule is not a certificate of global optimality, and a measured NC/TC accuracy difference must not be presented as a theorem of intrinsic statistical inferiority.\n\n');
fprintf(fid,'The stationary/seasonal comparison uses common error realizations when signal-error correlation is zero. With unit gains, truth cancels in pairwise differences; nearly identical TCH fitted variance results can therefore be structural. Its first-output full-covariance reconstruction can still vary with Q, and severe near-collinearity can amplify optimizer roundoff sensitivity.\n\n');
fprintf(fid,'All results describe the supplied synthetic generator. Do not claim universally superior accuracy, universal failure-free behavior, identified physical error covariances, or completed global validation from these outputs.\n');
end

function s=fieldnum(st,name)
if ~isfield(st,name), s='NA'; return; end
v=st.(name);
if isscalar(v) && isfinite(v) && v==fix(v), s=sprintf('%.0f',v); else, s=num(v); end
end

function write_guide(path,files,cfg)
fid=fopen(path,'w'); assert(fid>=0); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Synthetic experiment artifacts\n\n');
fprintf(fid,'Generated by `Build_Synthetic_Revision_Artifacts` from existing numeric outputs. This step does not fit models, select favorable results, or execute global experiments. PNG files are 600 dpi and matching PDF files are vector exports.\n\n');
fprintf(fid,'- Tables R1-R4: experiment design, complete scenario performance, paired temporal changes, and common-valid-subset method contrasts.\n');
fprintf(fid,'- Table R5 and Figure R3: severe-correlation stress, finite-covariance denominators, native PSD correction, and the paired test-MSE change relative to raw NC with the same downstream safeguard.\n');
fprintf(fid,'- Table R6 and Figures R4/R4b: eigenvalue-floor and initialization sensitivity, respectively.\n');
fprintf(fid,'- Tables R7-R8: p=2 nonuniqueness and p=3 physical-selection ambiguity.\n');
fprintf(fid,'- Figure R1 files: variance estimation at each configured signal-error correlation, using non-stress base correlations when available. The main plots show TC, TCH and NC PSD; full raw-NC diagnostics remain in the tables and stress figure.\n');
fprintf(fid,'- Figure R2: held-out merging performance for TC, TCH and NC PSD. Figure R5 shows NC PSD and the two controls; Figure R6 shows paired seasonal changes. All estimated-covariance methods use the same WA safeguard. Bars show the source 95%% Monte Carlo mean CIs; unavailable CIs have no drawn bar and remain NA in tables/text.\n\n');
fprintf(fid,'NA means unavailable or no valid denominator. It does not mean zero error, zero failure risk, or identical method performance. If only one sensitivity setting or no relevant scenario was supplied, the affected panel is explicitly unavailable.\n\n');
if isfield(cfg,'ciMethod'), fprintf(fid,'Source CI rule: %s.\n\n',char(cfg.ciMethod)); end
fprintf(fid,'## Files\n\n');
for k=1:numel(files), [~,name,ext]=fileparts(files{k}); fprintf(fid,'- `%s%s`\n',name,ext); end
fprintf(fid,'\nFigure/table identifiers are code identifiers. See `synthetic_results.md` for automatically populated counts, contrasts and limitations, and the repository README for figure mappings.\n');
end
