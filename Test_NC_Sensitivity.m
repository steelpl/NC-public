clear; clc
%% ============================================================
% Test_NC_Sensitivity.m
%
% Additional sensitivity tests for N-Tuple Collocation (NC)
%
% Experiments:
%   (1) SNR sensitivity:
%       SNRdB = [0.1 1 10]
%       p = 3
%       ecc = [0.3 0.6]
%       a = random, normalized to mean one
%       input covariance = true population covariance Q
%
%   (2) Number-of-datasets sensitivity:
%       p = 2:2:10
%       SNRdB = 0.1
%       ecc = [0.3 0.6]
%       a = random, normalized to mean one
%       input covariance = true population covariance Q
%
%   (3) Finite-sample sensitivity:
%       n = [10 30 50 100 200]
%       p = 3
%       SNRdB = 0.1
%       ecc = [0.3 0.6]
%       a = random, normalized to mean one
%       input covariance = sample covariance Qhat estimated from X
%
% Metric:
%   - MSE of estimated error variance
%
% Figure is saved to:
%   figures/Test_NC_Sensitivity_results.png
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
%% ============================================================

%% ---------------- Parameters ----------------
rng(123);

nMC = 500;
tol = 1e-10;

mean_a = 1.0;
ecc_list = [0.3 0.6];

SNRdB_list = [0.1 1 10];
p_fixed = 3;

p_list = 2:2:10;
SNRdB_fixed = 0.1;

n_list = [10 30 50 100 200];

%% ---------------- Run experiments ----------------
results_snr = run_snr_sensitivity( ...
    p_fixed, nMC, SNRdB_list, ecc_list, mean_a, tol);

results_p = run_p_sensitivity( ...
    p_list, nMC, SNRdB_fixed, ecc_list, mean_a, tol);

results_n = run_n_sensitivity( ...
    p_fixed, nMC, n_list, SNRdB_fixed, ecc_list, mean_a, tol);

%% ---------------- Plot ----------------
plot_NCtest_results(results_snr, results_p, results_n, ...
    SNRdB_list, p_list, n_list, ecc_list)


%% ============================================================
% Local function 1: SNR sensitivity
%% ============================================================
function R = run_snr_sensitivity(p,nMC,SNRdB_list,ecc_list,mean_a,tol)

nS = numel(SNRdB_list);
nE = numel(ecc_list);

mse_nc  = nan(nMC,nS,nE);
fail_nc = false(nMC,nS,nE);

% Common random numbers across SNR and ecc levels
A_all    = rand(p,nMC);
varE_all = rand(p,nMC);

for ie = 1:nE
    
    ecc = ecc_list(ie);
    
    for is = 1:nS
        
        SNRdB = SNRdB_list(is);
        SNR = 10^(SNRdB/10);
        
        for k = 1:nMC
            
            a = generate_scaling(A_all(:,k),mean_a);
            
            [E_true,Q] = generate_covariance_system( ...
                p,ecc,varE_all(:,k),SNR,a);
            
            [E_NC,~,~] = NC(Q);
            
            fail_nc(k,is,ie) = is_invalid_EeeT(E_NC,tol);
            
            if ~fail_nc(k,is,ie)
                mse_nc(k,is,ie) = error_variance_mse(E_NC,E_true);
            end
            
        end
        
    end
    
end

R = summarize_experiment(mse_nc,fail_nc);

end


%% ============================================================
% Local function 2: p sensitivity
%% ============================================================
function R = run_p_sensitivity(p_list,nMC,SNRdB,ecc_list,mean_a,tol)

nP = numel(p_list);
nE = numel(ecc_list);

SNR = 10^(SNRdB/10);

mse_nc  = nan(nMC,nP,nE);
fail_nc = false(nMC,nP,nE);

for ie = 1:nE
    
    ecc = ecc_list(ie);
    
    for ip = 1:nP
        
        p = p_list(ip);
        
        % Common random numbers within each p and ecc case
        A_all    = rand(p,nMC);
        varE_all = rand(p,nMC);
        
        for k = 1:nMC
            
            a = generate_scaling(A_all(:,k),mean_a);
            
            [E_true,Q] = generate_covariance_system( ...
                p,ecc,varE_all(:,k),SNR,a);
            
            [E_NC,~,~] = NC(Q);
            
            fail_nc(k,ip,ie) = is_invalid_EeeT(E_NC,tol);
            
            if ~fail_nc(k,ip,ie)
                mse_nc(k,ip,ie) = error_variance_mse(E_NC,E_true);
            end
            
        end
        
    end
    
end

R = summarize_experiment(mse_nc,fail_nc);

end


%% ============================================================
% Local function 3: finite-sample sensitivity
%% ============================================================
function R = run_n_sensitivity(p,nMC,n_list,SNRdB,ecc_list,mean_a,tol)

nN = numel(n_list);
nE = numel(ecc_list);

SNR = 10^(SNRdB/10);

mse_nc  = nan(nMC,nN,nE);
fail_nc = false(nMC,nN,nE);

% Common random numbers for true covariance generation
A_all    = rand(p,nMC);
varE_all = rand(p,nMC);

for ie = 1:nE
    
    ecc = ecc_list(ie);
    
    for in = 1:nN
        
        n = n_list(in);
        
        for k = 1:nMC
            
            a = generate_scaling(A_all(:,k),mean_a);
            
            [E_true,Q] = generate_covariance_system( ...
                p,ecc,varE_all(:,k),SNR,a);
            
            % Generate finite samples and estimate sample covariance Qhat.
            X = generate_centered_samples_from_cov(Q,n);
            Qhat = (X' * X) / (n - 1);
            Qhat = 0.5 * (Qhat + Qhat');   % suppress numerical asymmetry
            
            [E_NC,~,~] = NC(Qhat);
            
            fail_nc(k,in,ie) = is_invalid_EeeT(E_NC,tol);
            
            if ~fail_nc(k,in,ie)
                mse_nc(k,in,ie) = error_variance_mse(E_NC,E_true);
            end
            
        end
        
    end
    
end

R = summarize_experiment(mse_nc,fail_nc);

end


%% ============================================================
% Local function 4: Generate covariance system
%% ============================================================
function [E_true,Q] = generate_covariance_system(p,ecc,var_e,SNR,a)

EeeT_true = EeeTGEN(p,ecc,var_e);
E_true = diag(EeeT_true);

sigma_y2 = mean(E_true*SNR);
Q = EeeT_true + sigma_y2*(a*a');
Q = 0.5 * (Q + Q');

end


%% ============================================================
% Local function 5: Generate centered samples from covariance
%% ============================================================
function X = generate_centered_samples_from_cov(Q,n)

p = size(Q,1);

% Symmetrize and stabilize Q.
Q = 0.5 * (Q + Q');

% Eigen-based square root avoids mvnrnd and reduces toolbox dependency.
[V,D] = eig(Q);
d = diag(D);
d(d < 0 & d > -1e-12) = 0;

if any(d < 0)
    error('Q is not positive semidefinite beyond numerical tolerance.');
end

A = V * diag(sqrt(d));

Z = randn(n,p);
X = Z * A';

% Center samples before covariance calculation.
X = X - mean(X,1);

end


%% ============================================================
% Local function 6: Generate scaling vector
%% ============================================================
function a = generate_scaling(a_raw,mean_a)

a = 2*a_raw;
a = a/mean(a);     % normalize to mean(a)=1
a = mean_a*a;      % enforce target mean(a)

end


%% ============================================================
% Local function 7: Error variance MSE
%% ============================================================
function val = error_variance_mse(E_est,E_true)

val = mean((diag(E_est)-E_true).^2);

end


%% ============================================================
% Local function 8: Validity check
%% ============================================================
function tf = is_invalid_EeeT(E,tol)

tf = any(isnan(E(:))) || ...
     any(isinf(E(:))) || ...
     any(diag(E)<-tol);

end


%% ============================================================
% Local function 9: Summarize experiment
%% ============================================================
function R = summarize_experiment(mse_nc,fail_nc)

nX = size(mse_nc,2);
nE = size(mse_nc,3);

R.mse_med   = nan(nX,nE);
R.mse_q1    = nan(nX,nE);
R.mse_q3    = nan(nX,nE);
R.fail_rate = nan(nX,nE);

for ie = 1:nE
    [R.mse_med(:,ie),R.mse_q1(:,ie),R.mse_q3(:,ie)] = ...
        summarize_matrix(mse_nc(:,:,ie));
    
    R.fail_rate(:,ie) = squeeze(mean(fail_nc(:,:,ie),1)) * 100;
end

end


%% ============================================================
% Local function 10: Summary statistics for matrix columns
%% ============================================================
function [med,q1,q3] = summarize_matrix(X)

n = size(X,2);

med = nan(n,1);
q1  = nan(n,1);
q3  = nan(n,1);

for i = 1:n
    x = X(:,i);
    x = x(~isnan(x));
    [med(i),q1(i),q3(i)] = summarize_vector(x);
end

end


%% ============================================================
% Local function 11: Summary statistics for vector
%% ============================================================
function [med,q1,q3] = summarize_vector(x)

if isempty(x)
    med = nan;
    q1  = nan;
    q3  = nan;
else
    med = median(x);
    q1  = prctile(x,25);
    q3  = prctile(x,75);
end

end


%% ============================================================
% Local function 12: Manuscript-ready plot
%% ============================================================
function plot_NCtest_results(results_snr,results_p,results_n, ...
    SNRdB_list,p_list,n_list,ecc_list)

fig = figure('Units','centimeters', ...
             'Position',[2 2 24 7.8], ...
             'Color','w', ...
             'Renderer','painters');

set(fig,'PaperPositionMode','auto')

t = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

colors_fill = [
    0.70 0.80 1.00   % ecc = 0.3
    1.00 0.75 0.70   % ecc = 0.6
];

colors_line = [
    0.00 0.00 1.00   % ecc = 0.3
    1.00 0.00 0.00   % ecc = 0.6
];

legend_labels = compose('ecc = %.1f',ecc_list);

%% ---------- (a) SNR sensitivity ----------
ax = nexttile(t,1);
hold(ax,'on')

x_snr = 1:numel(SNRdB_list);

for ie = 1:numel(ecc_list)
    plot_iqr_shading(ax,x_snr, ...
        results_snr.mse_med(:,ie), ...
        results_snr.mse_q1(:,ie), ...
        results_snr.mse_q3(:,ie), ...
        colors_fill(ie,:),colors_line(ie,:))
end

xlabel(ax,'SNR_{dB}')
ylabel(ax,'Error Variance MSE')
title(ax,'(a) SNR sensitivity | p = 3')

xticks(ax,x_snr)
xticklabels(ax,compose('%.1g',SNRdB_list))
xlim(ax,[0.75 numel(SNRdB_list)+0.25])

format_axes(ax)

%% ---------- (b) p sensitivity ----------
ax = nexttile(t,2);
hold(ax,'on')

for ie = 1:numel(ecc_list)
    plot_iqr_shading(ax,p_list, ...
        results_p.mse_med(:,ie), ...
        results_p.mse_q1(:,ie), ...
        results_p.mse_q3(:,ie), ...
        colors_fill(ie,:),colors_line(ie,:))
end

xlabel(ax,'Number of datasets, p')
ylabel(ax,'Error Variance MSE')
title(ax,'(b) Dataset-number sensitivity | SNR_{dB} = 0.1')

xticks(ax,p_list)
xlim(ax,[min(p_list)-0.5 max(p_list)+0.8])

format_axes(ax)

%% ---------- (c) finite-sample sensitivity ----------
ax = nexttile(t,3);
hold(ax,'on')

x_n = 1:numel(n_list);

h_lines = gobjects(numel(ecc_list),1);

for ie = 1:numel(ecc_list)
    h_lines(ie) = plot_iqr_shading(ax,x_n, ...
        results_n.mse_med(:,ie), ...
        results_n.mse_q1(:,ie), ...
        results_n.mse_q3(:,ie), ...
        colors_fill(ie,:),colors_line(ie,:));
end

xlabel(ax,'Number of samples, n')
ylabel(ax,'Error Variance MSE')
title(ax,'(c) Finite-sample sensitivity | p = 3')

xticks(ax,x_n)
xticklabels(ax,string(n_list))
xlim(ax,[0.75 numel(n_list)+0.25])

format_axes(ax)

legend(ax,h_lines,legend_labels, ...
    'Location','northeast', ...
    'Box','on')

%% ---------- Harmonize y-axis limits ----------
ymax = max([
    results_snr.mse_q3(:);
    results_p.mse_q3(:);
    results_n.mse_q3(:)
],[],'omitnan');

ymax = ceil(ymax*10)/10;
if ymax <= 0
    ymax = 0.1;
end

axs = findall(fig,'Type','axes');
for i = 1:numel(axs)
    ylim(axs(i),[0 ymax])
end

%% ---------- Save figure ----------
save_figure(fig,'Test_NC_Sensitivity_results.png')

end


%% ============================================================
% Local function 13: Plot median + IQR
%% ============================================================
function h = plot_iqr_shading(ax,x,med,q1,q3,shade_color,line_color)

hold(ax,'on')

fill(ax, ...
    [x(:)' fliplr(x(:)')], ...
    [q1(:)' fliplr(q3(:)')], ...
    shade_color, ...
    'EdgeColor','none', ...
    'FaceAlpha',0.30, ...
    'HandleVisibility','off')

h = plot(ax,x(:),med(:), ...
    '-o', ...
    'Color',line_color, ...
    'MarkerFaceColor',line_color, ...
    'MarkerEdgeColor',line_color, ...
    'MarkerSize',4.2, ...
    'LineWidth',1.7);

end


%% ============================================================
% Local function 14: Axis formatting
%% ============================================================
function format_axes(ax)

grid(ax,'on')
box(ax,'on')

set(ax, ...
    'FontName','Arial', ...
    'FontSize',9.0, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Layer','top')

ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

end


%% ============================================================
% Local function 15: Save figure
%% ============================================================
function save_figure(fig,filename)

outDir = fullfile(pwd,'figures');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

exportgraphics(fig,fullfile(outDir,filename), ...
    'Resolution',600, ...
    'BackgroundColor','white');

[~,name,~] = fileparts(filename);
exportgraphics(fig,fullfile(outDir,[name,'.pdf']), ...
    'ContentType','vector', ...
    'BackgroundColor','white');

end
