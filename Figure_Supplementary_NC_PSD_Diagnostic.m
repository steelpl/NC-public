clear; close all; clc;

%% ============================================================
% Figure_Supplementary_NC_PSD_Diagnostic.m
%
% Purpose:
%   Supplementary robustness check for optional PSD-projected NC.
%
% Diagnostics:
%   1) Fraction of synthetic pixels requiring PSD correction
%   2) Theta-shrinkage factor distribution
%   3) Objective increase after PSD correction
%   4) WA validation-proxy ubRMSE/R before vs after PSD correction
%
% Notes:
%   NC defaults to PSD enforcement. This script explicitly compares raw NC
%   with its PSD feasibility map as a supplementary numerical diagnostic.
%   Its validation proxy is evaluated on the same synthetic samples used
%   for estimation; use run_revision_synthetic for independent-test metrics.
%% ============================================================

%% User settings

rng(20260707);

nPixels = 1000;
nObs = 365;
pList = [3 4 5 6];
minEigTol = 1e-10;

outRoot = fullfile(pwd, 'figures', 'Supplementary_NC_PSD_Diagnostic');
figFld = fullfile(outRoot, 'figures');
tableFld = fullfile(outRoot, 'tables');
makeFolder(outRoot);
makeFolder(figFld);
makeFolder(tableFld);

exportPNG = true;
exportPDF = true;
figureDPI = 600;
fontName = 'Arial';
fontSize = 12;

%% Run diagnostic simulation

D = runPSDProjectionDiagnostic(nPixels, nObs, pList, minEigTol);
S = summarizeDiagnostics(D);

writetable(D, fullfile(tableFld, 'NC_PSD_Diagnostic_PixelLevel.csv'));
writetable(S, fullfile(tableFld, 'NC_PSD_Diagnostic_Summary.csv'));

fprintf('============================================================\n');
fprintf('Supplementary NC PSD diagnostic\n');
fprintf('Synthetic pixels: %d\n', height(D));
fprintf('PSD correction required: %.2f %%\n', 100 * S.Value(strcmp(S.Metric, 'fractionRequiringPSD'))); 
fprintf('Median theta shrinkage alpha, all/corrected: %.4f / %.4f\n', ...
    S.Value(strcmp(S.Metric, 'medianThetaScale')), S.Value(strcmp(S.Metric, 'medianThetaScale_corrected'))); 
fprintf('Median relative objective increase, all/corrected: %.4f / %.4f %%\n', ...
    100 * S.Value(strcmp(S.Metric, 'medianRelativeObjectiveIncrease')), ...
    100 * S.Value(strcmp(S.Metric, 'medianRelativeObjectiveIncrease_corrected'))); 
fprintf('Median ubRMSE raw / PSD: %.4f / %.4f\n', ...
    S.Value(strcmp(S.Metric, 'medianUbRMSE_raw')), S.Value(strcmp(S.Metric, 'medianUbRMSE_psd')));
fprintf('Median R raw / PSD: %.4f / %.4f\n', ...
    S.Value(strcmp(S.Metric, 'medianR_raw')), S.Value(strcmp(S.Metric, 'medianR_psd')));
fprintf('Output folder: %s\n', outRoot);
fprintf('============================================================\n');

%% Plot core supplementary figure

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 19 14]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plotCorrectionFraction(nexttile(tl, 1), D, fontName, fontSize);
plotThetaShrinkage(nexttile(tl, 2), D, fontName, fontSize);
plotObjectiveIncrease(nexttile(tl, 3), D, fontName, fontSize);
plotValidationBeforeAfter(nexttile(tl, 4), D, fontName, fontSize);

exportFigure(fig, figFld, 'FigS_NC_PSD_Diagnostic', figureDPI, exportPNG, exportPDF);

fprintf('Finished. Figure saved in:\n%s\n', figFld);

%% Local functions

function D = runPSDProjectionDiagnostic(nPixels, nObs, pList, minEigTol)

rows = repmat(emptyRow(), nPixels, 1);

optsRaw = struct('enforcePSD', false, 'repairInputPSD', true, ...
    'nStarts', 4, 'nRandomStarts', 0, 'verbose', false);
optsPSD = optsRaw;
optsPSD.enforcePSD = true;

for k = 1:nPixels
    p = pList(randi(numel(pList)));
    ecc = 0.20 + 0.75 * rand;
    snrDB = -5 + 10 * rand;
    snr = 10^(snrDB / 10);
    meanA = 0.75 + 0.50 * rand;

    a = randomScalingVector(p, meanA);
    varE = 0.05 + 0.95 * rand(p, 1);
    Etrue = EeeTGEN(p, ecc, varE);
    signalVar = mean(diag(Etrue) * snr);

    [X, signal] = simulateObservationSystem(nObs, a, signalVar, Etrue);
    Q = cov(X, 1);

    [Eraw, ~, ~, infoRaw] = NC(Q, optsRaw);
    [Epsd, ~, ~, infoPSD] = NC(Q, optsPSD);

    minEigRaw = min(eig(0.5 * (Eraw + Eraw')));
    minEigPSD = min(eig(0.5 * (Epsd + Epsd')));
    scaleRaw = max(1, norm(Eraw, 'fro'));
    requiresPSD = minEigRaw < -minEigTol * scaleRaw;

    [ubRaw, rRaw, validRaw] = evaluateWAValidationProxy(Eraw, X, signal);
    [ubPSD, rPSD, validPSD] = evaluateWAValidationProxy(Epsd, X, signal);

    rows(k).pixelID = k;
    rows(k).p = p;
    rows(k).ecc = ecc;
    rows(k).snrDB = snrDB;
    rows(k).meanA = meanA;
    rows(k).requiresPSD = requiresPSD;
    rows(k).thetaScalePSD = infoPSD.thetaScalePSD;
    rows(k).objectiveBeforePSD = infoPSD.objectiveBeforePSD;
    rows(k).objectiveAfterPSD = infoPSD.objective;
    rows(k).relativeObjectiveIncrease = max(0, infoPSD.objective - infoPSD.objectiveBeforePSD) ./ ...
        max(infoPSD.objectiveBeforePSD, eps);
    rows(k).minEigRaw = minEigRaw;
    rows(k).minEigPSD = minEigPSD;
    rows(k).rawFeasiblePSD = infoRaw.feasiblePSD;
    rows(k).psdFeasiblePSD = infoPSD.feasiblePSD;
    rows(k).ubRMSE_raw = ubRaw;
    rows(k).ubRMSE_psd = ubPSD;
    rows(k).R_raw = rRaw;
    rows(k).R_psd = rPSD;
    rows(k).validWA_raw = validRaw;
    rows(k).validWA_psd = validPSD;
end

D = struct2table(rows);

end

function row = emptyRow()
row = struct('pixelID', NaN, 'p', NaN, 'ecc', NaN, 'snrDB', NaN, 'meanA', NaN, ...
    'requiresPSD', false, 'thetaScalePSD', NaN, ...
    'objectiveBeforePSD', NaN, 'objectiveAfterPSD', NaN, ...
    'relativeObjectiveIncrease', NaN, 'minEigRaw', NaN, 'minEigPSD', NaN, ...
    'rawFeasiblePSD', false, 'psdFeasiblePSD', false, ...
    'ubRMSE_raw', NaN, 'ubRMSE_psd', NaN, 'R_raw', NaN, 'R_psd', NaN, ...
    'validWA_raw', false, 'validWA_psd', false);
end

function a = randomScalingVector(p, meanA)
a = 0.4 + 1.2 * rand(p, 1);
a = a / mean(a) * meanA;
end

function [X, signal] = simulateObservationSystem(nObs, a, signalVar, Etrue)
signal = sqrt(signalVar) * randn(nObs, 1);
Ehalf = psdSquareRoot(Etrue);
err = randn(nObs, numel(a)) * Ehalf';
X = signal * a' + err;
end

function Ahalf = psdSquareRoot(A)
A = 0.5 * (A + A');
[V, D] = eig(A);
d = real(diag(D));
d(d < 0) = 0;
Ahalf = V * diag(sqrt(d));
end

function [ubrmse, rValue, isValid] = evaluateWAValidationProxy(Eest, X, reference)
ubrmse = NaN;
rValue = NaN;
isValid = false;
try
    w = WA(Eest);
    merged = X * w;
    ok = isfinite(merged) & isfinite(reference);
    if nnz(ok) < 30
        return;
    end
    merged = merged(ok);
    reference = reference(ok);
    d = merged - reference;
    d = d - mean(d, 'omitnan');
    ubrmse = sqrt(mean(d.^2, 'omitnan'));
    C = corrcoef(merged, reference, 'Rows', 'complete');
    if numel(C) >= 4
        rValue = C(1, 2);
    end
    isValid = isfinite(ubrmse) && isfinite(rValue);
catch
    isValid = false;
end
end

function S = summarizeDiagnostics(D)
corrected = D.requiresPSD;
metrics = [
    "fractionRequiringPSD"
    "medianThetaScale"
    "medianThetaScale_corrected"
    "medianRelativeObjectiveIncrease"
    "medianRelativeObjectiveIncrease_corrected"
    "medianUbRMSE_raw"
    "medianUbRMSE_psd"
    "medianR_raw"
    "medianR_psd"
    "fractionWAValid_raw"
    "fractionWAValid_psd"
    ];
values = [
    mean(D.requiresPSD)
    median(D.thetaScalePSD, 'omitnan')
    medianOrNaN(D.thetaScalePSD(corrected))
    median(D.relativeObjectiveIncrease, 'omitnan')
    medianOrNaN(D.relativeObjectiveIncrease(corrected))
    median(D.ubRMSE_raw, 'omitnan')
    median(D.ubRMSE_psd, 'omitnan')
    median(D.R_raw, 'omitnan')
    median(D.R_psd, 'omitnan')
    mean(D.validWA_raw)
    mean(D.validWA_psd)
    ];
S = table(metrics, values, 'VariableNames', {'Metric', 'Value'});
end

function value = medianOrNaN(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function plotCorrectionFraction(ax, D, fontName, fontSize)
frac = mean(D.requiresPSD);
b = bar(ax, categorical({'PSD correction'}), frac, 0.55);
b.FaceColor = [0.25 0.42 0.70];
ylim(ax, [0 1]);
ylabel(ax, 'Fraction of synthetic pixels');
title(ax, '(a) PSD correction frequency');
grid(ax, 'on');
text(ax, 1, min(frac + 0.04, 0.96), sprintf('%.1f%%', 100 * frac), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
formatAxes(ax, fontName, fontSize);
end

function plotThetaShrinkage(ax, D, fontName, fontSize)
alpha = D.thetaScalePSD(D.requiresPSD & isfinite(D.thetaScalePSD));
if isempty(alpha)
    text(ax, 0.5, 0.5, 'No PSD correction required', 'HorizontalAlignment', 'center');
    xlim(ax, [0 1]);
else
    histogram(ax, alpha, min(16, max(4, numel(unique(alpha)))), ...
        'FaceColor', [0.18 0.58 0.42], 'EdgeColor', 'none');
    xlim(ax, [max(0, min(alpha) - 0.02) 1]);
    text(ax, 0.03, 0.92, sprintf('Corrected pixels: n = %d', numel(alpha)), ...
        'Units', 'normalized', 'FontWeight', 'bold');
end
xlabel(ax, 'Theta-shrinkage factor alpha');
ylabel(ax, 'Count');
title(ax, '(b) Correction strength');
grid(ax, 'on');
formatAxes(ax, fontName, fontSize);
end

function plotObjectiveIncrease(ax, D, fontName, fontSize)
value = 100 * D.relativeObjectiveIncrease(D.requiresPSD & isfinite(D.relativeObjectiveIncrease));
value = value(value >= 0);
if isempty(value)
    text(ax, 0.5, 0.5, 'No PSD correction required', 'HorizontalAlignment', 'center');
else
    histogram(ax, value, min(18, max(4, numel(unique(round(value, 6))))), ...
        'FaceColor', [0.78 0.43 0.22], 'EdgeColor', 'none');
    text(ax, 0.03, 0.92, sprintf('Median = %.2f%%', median(value, 'omitnan')), ...
        'Units', 'normalized', 'FontWeight', 'bold');
end
xlabel(ax, 'Objective increase after PSD correction (%)');
ylabel(ax, 'Count');
title(ax, '(c) NC fit degradation');
grid(ax, 'on');
formatAxes(ax, fontName, fontSize);
end

function plotValidationBeforeAfter(ax, D, fontName, fontSize)
valid = D.validWA_raw & D.validWA_psd;

ubRaw = median(D.ubRMSE_raw(valid), 'omitnan');
ubPSD = median(D.ubRMSE_psd(valid), 'omitnan');
rRaw = median(D.R_raw(valid), 'omitnan');
rPSD = median(D.R_psd(valid), 'omitnan');

hold(ax, 'on');
yyaxis(ax, 'left');
b1 = bar(ax, [1 2], [ubRaw ubPSD], 0.55);
b1.FaceColor = [0.35 0.35 0.35];
ylabel(ax, 'Median ubRMSE');
ylim(ax, [0 max([ubRaw ubPSD]) * 1.35 + eps]);

yyaxis(ax, 'right');
b2 = bar(ax, [4 5], [rRaw rPSD], 0.55);
b2.FaceColor = [0.25 0.55 0.78];
ylabel(ax, 'Median R');
ylim(ax, [0 1]);

xlim(ax, [0.3 5.7]);
xticks(ax, [1 2 4 5]);
xticklabels(ax, {'ubRMSE raw', 'ubRMSE PSD', 'R raw', 'R PSD'});
xtickangle(ax, 25);
title(ax, '(d) Validation proxy before vs after');
grid(ax, 'on');
formatAxes(ax, fontName, fontSize);
end

function formatAxes(ax, fontName, fontSize)
set(ax, 'FontName', fontName, 'FontSize', fontSize, 'Box', 'off', ...
    'LineWidth', 0.8, 'TickDir', 'out');
end

function makeFolder(folderPath)
if exist(folderPath, 'dir') ~= 7
    mkdir(folderPath);
end
end

function exportFigure(fig, figFld, baseName, dpi, exportPNG, exportPDF)
if exportPNG
    exportgraphics(fig, fullfile(figFld, [baseName, '.png']), 'Resolution', dpi);
end
if exportPDF
    exportgraphics(fig, fullfile(figFld, [baseName, '.pdf']), 'ContentType', 'vector');
end
end
