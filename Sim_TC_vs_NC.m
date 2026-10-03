clear; clc
%% ============================================================
% TC vs NC analysis
%
% Outputs:
%   (1) Difference (TC - NC) under increasing ecc
%   (2) MSE of estimated error variance under increasing ecc
%   (3) Failure rate under increasing ecc
%   (4) Scaling sensitivity under varying mean(a)
%
% Figure is saved to:
%   figures/Sim_TC_vs_NC_results.png
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
%% ============================================================

%% ---------------- Parameters ----------------
p = 3;
rng(123);

ecc_list    = [0:0.1:0.90, 0.99];  % error cross-correlation levels
mean_a_list = 1:0.5:5;              % scaling sensitivity levels
ecc_fixed   = 0.3;                  % fixed ecc for scaling test

SNRdB = 0.1;
SNR   = 10^(SNRdB/10);

nMC = 500;                          % Monte Carlo runs
tol = 1e-10;                        % numerical tolerance

%% ---------------- Common random inputs ----------------
A_all    = rand(p,nMC);             % scaling seeds
varE_all = rand(p,nMC);             % error variance seeds

%% ---------------- Simulations ----------------
results_ecc = run_TCvsNC_ecc_sweep( ...
    p, nMC, ecc_list, SNR, A_all, varE_all, tol);

results_scale = run_TCvsNC_scaling_sensitivity( ...
    p, nMC, ecc_fixed, mean_a_list, SNR, A_all, varE_all, tol);

%% ---------------- Plot ----------------
plot_TCvsNC_results(results_ecc, results_scale, ecc_list, mean_a_list, ecc_fixed);


%% ============================================================
% Local function 1: ECC sweep
%% ============================================================
function R = run_TCvsNC_ecc_sweep(p,nMC,ecc_list,SNR,A_all,varE_all,tol)

nE = numel(ecc_list);

diff_all = nan(nMC,p,nE);
mse_tc   = nan(nMC,nE);
mse_nc   = nan(nMC,nE);

fail_tc = false(nMC,nE);
fail_nc = false(nMC,nE);

for ie = 1:nE
    
    ecc = ecc_list(ie);
    
    for k = 1:nMC
        
        a = generate_scaling(A_all(:,k),1);   % random a with mean(a)=1
        
        [E_true,ExxT] = generate_covariance_system( ...
            p,ecc,varE_all(:,k),SNR,a);
        
        [E_TC,E_NC] = estimate_TC_NC(ExxT);
        
        fail_tc(k,ie) = is_invalid_EeeT(E_TC,tol);
        fail_nc(k,ie) = is_invalid_EeeT(E_NC,tol);
        
        if ~fail_tc(k,ie)
            mse_tc(k,ie) = error_variance_mse(E_TC,E_true);
        end
        
        if ~fail_nc(k,ie)
            mse_nc(k,ie) = error_variance_mse(E_NC,E_true);
        end
        
        if ~fail_tc(k,ie) && ~fail_nc(k,ie)
            diff_all(k,:,ie) = diag(E_TC) - diag(E_NC);
        end
        
    end
    
end

% Difference statistics across MC samples and diagonal elements
R.diff_med = nan(nE,1);
R.diff_q1  = nan(nE,1);
R.diff_q3  = nan(nE,1);

for ie = 1:nE
    x = diff_all(:,:,ie);
    x = x(~isnan(x));
    [R.diff_med(ie),R.diff_q1(ie),R.diff_q3(ie)] = summarize_vector(x);
end

% MSE statistics
[R.mse_tc_med,R.mse_tc_q1,R.mse_tc_q3] = summarize_matrix(mse_tc);
[R.mse_nc_med,R.mse_nc_q1,R.mse_nc_q3] = summarize_matrix(mse_nc);

% Failure rate
R.tc_fail_rate = mean(fail_tc,1)'*100;
R.nc_fail_rate = mean(fail_nc,1)'*100;

end


%% ============================================================
% Local function 2: Scaling sensitivity
%% ============================================================
function R = run_TCvsNC_scaling_sensitivity( ...
    p,nMC,ecc_fixed,mean_a_list,SNR,A_all,varE_all,tol)

nA = numel(mean_a_list);

mse_tc = nan(nMC,nA);
mse_nc = nan(nMC,nA);

for ia = 1:nA
    
    mean_a = mean_a_list(ia);
    
    for k = 1:nMC
        
        a = generate_scaling(A_all(:,k),mean_a);  % controlled mean(a)
        
        [E_true,ExxT] = generate_covariance_system( ...
            p,ecc_fixed,varE_all(:,k),SNR,a);
        
        [E_TC,E_NC] = estimate_TC_NC(ExxT);
        
        if ~is_invalid_EeeT(E_TC,tol)
            mse_tc(k,ia) = error_variance_mse(E_TC,E_true);
        end
        
        if ~is_invalid_EeeT(E_NC,tol)
            mse_nc(k,ia) = error_variance_mse(E_NC,E_true);
        end
        
    end
    
end

[R.mse_tc_med,R.mse_tc_q1,R.mse_tc_q3] = summarize_matrix(mse_tc);
[R.mse_nc_med,R.mse_nc_q1,R.mse_nc_q3] = summarize_matrix(mse_nc);

end


%% ============================================================
% Local function 3: Covariance system
%% ============================================================
function [E_true,ExxT] = generate_covariance_system(p,ecc,var_e,SNR,a)

EeeT_true = EeeTGEN(p,ecc,var_e);
E_true    = diag(EeeT_true);

Ey2  = mean(E_true*SNR);          % signal variance
ExxT = EeeT_true + Ey2*(a*a');    % observation covariance

end


%% ============================================================
% Local function 4: Estimate TC and NC
%% ============================================================
function [E_TC,E_NC] = estimate_TC_NC(ExxT)

[E_TC,~,~] = TC(ExxT);
[E_NC,~,~] = NC(ExxT);

end


%% ============================================================
% Local function 5: Scaling vector
%% ============================================================
function a = generate_scaling(a_raw,mean_a)

a = 2*a_raw;
a = a/mean(a);       % normalize to mean(a)=1
a = mean_a*a;        % enforce target mean(a)

end


%% ============================================================
% Local function 6: Error variance MSE
%% ============================================================
function val = error_variance_mse(E_est,E_true)

val = mean((diag(E_est)-E_true).^2);

end


%% ============================================================
% Local function 7: Invalid EeeT check
%% ============================================================
function tf = is_invalid_EeeT(E,tol)

tf = any(isnan(E(:))) || ...
     any(isinf(E(:))) || ...
     any(diag(E)<-tol);

end


%% ============================================================
% Local function 8: Summary statistics for matrix columns
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
% Local function 9: Summary statistics for vector
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
% Local function 10: Manuscript-ready plot
%% ============================================================
function plot_TCvsNC_results(results_ecc,results_scale,ecc_list,mean_a_list,ecc_fixed)

fig = figure('Units','centimeters', ...
             'Position',[2 2 18 13], ...
             'Color','w');

t = tiledlayout(fig,2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Colors
col_tc_fill = [1.0 0.70 0.70];
col_tc_line = [1.0 0.00 0.00];
col_nc_fill = [0.70 0.80 1.00];
col_nc_line = [0.00 0.00 1.00];

%% ---------- (a) Difference ----------
ax = nexttile(t,1);

plot_iqr_shading(ax,ecc_list, ...
    results_ecc.diff_med*1e3, ...
    results_ecc.diff_q1*1e3, ...
    results_ecc.diff_q3*1e3, ...
    col_nc_fill,col_nc_line)

max_abs_diff = get_symmetric_limit([results_ecc.diff_q1; results_ecc.diff_q3]*1e3);

ylim(ax,[-max_abs_diff max_abs_diff])
yticks(ax,linspace(-max_abs_diff,max_abs_diff,5))

xlabel(ax,'ecc')
ylabel(ax,'TC - NC (\times10^{-3})')
title(ax,'(a) Difference (TC - NC) | random a')

format_axes(ax)
% add_panel_label(ax,'(a)')

%% ---------- (b) MSE under increasing ecc ----------
ax = nexttile(t,2);

plot_iqr_shading(ax,ecc_list, ...
    results_ecc.mse_tc_med, ...
    results_ecc.mse_tc_q1, ...
    results_ecc.mse_tc_q3, ...
    col_tc_fill,col_tc_line)

plot_iqr_shading(ax,ecc_list, ...
    results_ecc.mse_nc_med, ...
    results_ecc.mse_nc_q1, ...
    results_ecc.mse_nc_q3, ...
    col_nc_fill,col_nc_line)

ylim(ax,[0 0.5])
yticks(ax,0:0.1:0.5)

xlabel(ax,'ecc')
ylabel(ax,'Error Variance MSE')
title(ax,'(b) MSE | random a')

legend(ax, ...
    'TC IQR','TC median','NC IQR','NC median', ...
    'Location','northwest', ...
    'Box','on')

format_axes(ax)
% add_panel_label(ax,'(b)')

%% ---------- (c) Failure rate ----------
ax = nexttile(t,3);
hold(ax,'on')

plot(ax,ecc_list,results_ecc.tc_fail_rate, ...
    'o-', ...
    'Color',col_tc_line, ...
    'MarkerFaceColor',col_tc_line, ...
    'LineWidth',1.8, ...
    'MarkerSize',5)

plot(ax,ecc_list,results_ecc.nc_fail_rate, ...
    'o-', ...
    'Color',col_nc_line, ...
    'MarkerFaceColor',col_nc_line, ...
    'LineWidth',1.8, ...
    'MarkerSize',5)

xlabel(ax,'ecc')
ylabel(ax,'Failure rate (%)')
title(ax,'(c) Failure Rate | random a')

legend(ax,'TC','NC', ...
    'Location','northwest', ...
    'Box','on')

format_axes(ax)
% add_panel_label(ax,'(c)')

%% ---------- (d) Scaling sensitivity ----------
ax = nexttile(t,4);

plot_iqr_shading(ax,mean_a_list, ...
    results_scale.mse_tc_med, ...
    results_scale.mse_tc_q1, ...
    results_scale.mse_tc_q3, ...
    col_tc_fill,col_tc_line)

plot_iqr_shading(ax,mean_a_list, ...
    results_scale.mse_nc_med, ...
    results_scale.mse_nc_q1, ...
    results_scale.mse_nc_q3, ...
    col_nc_fill,col_nc_line)

ylim(ax,[0 0.3])
yticks(ax,0:0.1:0.3)

xlabel(ax,'mean(a)')
ylabel(ax,'Error Variance MSE')
title(ax,['(d) Scaling Sensitivity | ecc = ',num2str(ecc_fixed)])

% legend(ax, ...
%     'TC IQR','TC median','NC IQR','NC median', ...
%     'Location','northwest', ...
%     'Box','on')

format_axes(ax)
% add_panel_label(ax,'(d)')

save_figure(fig,'Sim_TC_vs_NC_results.png')

end


%% ============================================================
% Local function 11: Plot median + IQR
%% ============================================================
function plot_iqr_shading(ax,x,med,q1,q3,shade_color,line_color)

hold(ax,'on')

fill(ax, ...
    [x(:)' fliplr(x(:)')], ...
    [q1(:)' fliplr(q3(:)')], ...
    shade_color, ...
    'EdgeColor','none', ...
    'FaceAlpha',0.45)

plot(ax,x(:),med(:), ...
    '-', ...
    'Color',line_color, ...
    'LineWidth',2)

end


%% ============================================================
% Local function 12: Symmetric axis limit
%% ============================================================
function lim = get_symmetric_limit(x)

x = x(~isnan(x));

if isempty(x)
    lim = 1;
else
    lim = ceil(max(abs(x))*1.05);
    
    if lim == 0
        lim = 1;
    end
end

end


%% ============================================================
% Local function 13: Axis formatting
%% ============================================================
function format_axes(ax)

grid(ax,'on')
box(ax,'on')

set(ax, ...
    'FontName','Arial', ...
    'FontSize',10, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'Layer','top')

end


%% ============================================================
% Local function 14: Panel label
%% ============================================================
function add_panel_label(ax,labelText)

text(ax,0.02,0.96,labelText, ...
    'Units','normalized', ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top')

end


%% ============================================================
% Local function 15: Save figure
%% ============================================================
function save_figure(fig,filename)

outDir = fullfile(pwd,'figures');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

exportgraphics(fig,fullfile(outDir,filename),'Resolution',600);

end
