clear; clc
%% ============================================================
% WA merging performance using NC, TC, and TCH
%
% Purpose:
%   Compare WA fusion performance when the error covariance matrix is
%   estimated using:
%       1) NC
%       2) TC
%       3) TCH
%
% Metrics:
%   - MSE
%   - R^2
%
% Experiments:
%   1) ecc sweep:       ecc = [0.3 0.6 0.9], mean(a) = 1
%   2) mean(a) sweep:   mean(a) = [0.5 1.0 1.5], ecc = 0.3
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

%% ---------- Parameters ----------
p = 3;
rng(123);

SNRdB = 0.1;
SNR   = 10^(SNRdB/10);

nMC = 500;
tol = 1e-10;

method = {'NC','TC','TCH'};

method_colors = [
    0.00 0.60 0.00   % NC
    1.00 0.00 0.00   % TC
    0.00 0.00 1.00   % TCH
];

ecc_list = [0.3 0.6 0.9];
mean_a_fixed = 1.0;

ecc_fixed = 0.3;
mean_a_list = [0.5 1.0 1.5];

%% ---------- Common random inputs ----------
A_all    = rand(p,nMC);     % scaling seeds
varE_all = rand(p,nMC);     % error variance seeds

%% ---------- Simulations ----------
[MSE_ecc, R2_ecc] = run_comparison( ...
    p, nMC, ecc_list, mean_a_fixed, SNR, A_all, varE_all, tol, "ecc");

[MSE_a, R2_a] = run_comparison( ...
    p, nMC, mean_a_list, ecc_fixed, SNR, A_all, varE_all, tol, "mean_a");

%% ---------- Plot ----------
plot_WA_results( ...
    MSE_ecc, R2_ecc, MSE_a, R2_a, ...
    ecc_list, mean_a_list, ecc_fixed, method, method_colors);


%% ============================================================
% Local function 1: Core simulation
%% ============================================================
function [MSE_all,R2_all] = run_comparison( ...
    p,nMC,x_list,fixed_val,SNR,A_all,varE_all,tol,mode)

nX = numel(x_list);
nMethod = 3;   % NC, TC, TCH

MSE_all = nan(nMC,nX,nMethod);
R2_all  = nan(nMC,nX,nMethod);

for ix = 1:nX
    
    [ecc,mean_a] = resolve_case(x_list(ix),fixed_val,mode);
    
    for k = 1:nMC
        
        a = generate_scaling(A_all(:,k),mean_a);   % controlled scaling
        
        [ExxT,Ey2] = generate_covariance_system( ...
            p,ecc,varE_all(:,k),SNR,a);
        
        E_list = estimate_covariances(ExxT);       % order: NC, TC, TCH
        
        for m = 1:nMethod
            
            if is_valid_EeeT(E_list{m},tol)
                [MSE_all(k,ix,m),R2_all(k,ix,m)] = ...
                    evaluate_WA(E_list{m},ExxT,Ey2,a);
            end
            
        end
        
    end
end

end


%% ============================================================
% Local function 2: Resolve experiment case
%% ============================================================
function [ecc,mean_a] = resolve_case(x_value,fixed_val,mode)

switch mode
    case "ecc"
        ecc = x_value;
        mean_a = fixed_val;
        
    case "mean_a"
        ecc = fixed_val;
        mean_a = x_value;
        
    otherwise
        error('Unknown mode. Use "ecc" or "mean_a".')
end

end


%% ============================================================
% Local function 3: Scaling vector
%% ============================================================
function a = generate_scaling(a_raw,mean_a)

a = 2*a_raw;
a = a/mean(a);      % normalize to mean(a)=1
a = mean_a*a;       % enforce target mean

end


%% ============================================================
% Local function 4: Synthetic covariance system
%% ============================================================
function [ExxT,Ey2] = generate_covariance_system(p,ecc,var_e,SNR,a)

EeeT_true = EeeTGEN(p,ecc,var_e);

Ey2  = mean(diag(EeeT_true)*SNR);  % signal variance
ExxT = EeeT_true + Ey2*(a*a');     % observation covariance

end


%% ============================================================
% Local function 5: Estimate covariance matrices
%% ============================================================
function E_list = estimate_covariances(ExxT)

[E_NC,~,~]  = NC(ExxT);
[E_TC,~,~]  = TC(ExxT);
[E_TCH,~,~] = TCH(ExxT);

E_list = {E_NC,E_TC,E_TCH};   % method order

end


%% ============================================================
% Local function 6: WA performance
%% ============================================================
function [mse_val,r2_val] = evaluate_WA(EeeT,ExxT,Ey2,a)

try
    u = WA_Submitted(EeeT);
    mse_val = calc_MSE(u,ExxT,Ey2,a);
    r2_val  = calc_R2(u,ExxT,Ey2,a);
catch
    mse_val = nan;
    r2_val  = nan;
end

end


%% ============================================================
% Local function 7: MSE
%% ============================================================
function val = calc_MSE(u,ExxT,Ey2,a)

val = u'*ExxT*u - 2*Ey2*u'*a + Ey2;

end


%% ============================================================
% Local function 8: R2
%% ============================================================
function val = calc_R2(u,ExxT,Ey2,a)

den = max(u'*ExxT*u,1e-12);
val = Ey2*(u'*(a*a')*u) / den;

val = max(0,min(1,val));   % numerical cleanup

end


%% ============================================================
% Local function 9: Validity check
%% ============================================================
function tf = is_valid_EeeT(E,tol)

tf = true;

if any(isnan(E(:))) || any(isinf(E(:)))
    tf = false;
    return
end

if any(diag(E) < -tol)
    tf = false;
    return
end

E = 0.5*(E + E');

if min(eig(E)) < -tol
    tf = false;
end

end


%% ============================================================
% Local function 10: Manuscript-ready plot
%% ============================================================
function plot_WA_results( ...
    MSE_ecc,R2_ecc,MSE_a,R2_a, ...
    ecc_list,mean_a_list,ecc_fixed,method,method_colors)

fig = figure('Units','centimeters', ...
             'Position',[2 2 18 13], ...
             'Color','w');

t = tiledlayout(fig,2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% ---------- (a) MSE by ecc ----------
ax = nexttile(t,1);

plot_box_group(ax,MSE_ecc,ecc_list,method_colors,0.5,false)

ylabel(ax,'MSE')
xlabel(ax,'ecc')
title(ax,'(a) MSE by ecc | a = random')

format_axes(ax)
% add_panel_label(ax,'(a)')
add_method_legend(ax,method,method_colors,'northwest')

%% ---------- (b) R2 by ecc ----------
ax = nexttile(t,2);

plot_box_group(ax,R2_ecc,ecc_list,method_colors,0.2,true)

ylabel(ax,'R^2')
xlabel(ax,'ecc')
title(ax,'(b) R^2 by ecc | a = random')

format_axes(ax)
% add_panel_label(ax,'(b)')

%% ---------- (c) MSE by mean(a) ----------
ax = nexttile(t,3);

plot_box_group(ax,MSE_a,mean_a_list,method_colors,0.5,false)

ylabel(ax,'MSE')
xlabel(ax,'mean(a)')
title(ax,['(c) MSE by mean(a) | ecc = ',num2str(ecc_fixed)])

format_axes(ax)
% add_panel_label(ax,'(c)')

%% ---------- (d) R2 by mean(a) ----------
ax = nexttile(t,4);

plot_box_group(ax,R2_a,mean_a_list,method_colors,0.2,true)

ylabel(ax,'R^2')
xlabel(ax,'mean(a)')
title(ax,['(d) R^2 by mean(a) | ecc = ',num2str(ecc_fixed)])

format_axes(ax)
% add_panel_label(ax,'(d)')

save_figure(fig,'Sim_WA_Comparison.png')

end


%% ============================================================
% Local function 11: Grouped boxplot
%% ============================================================
function plot_box_group(ax,X,x_list,colors,ytick_step,isR2)

if nargin < 6
    isR2 = false;
end

axes(ax)
hold(ax,'on')

[~,nX,nMethod] = size(X);

box_width = 0.55;
box_pos = get_box_positions(nX,nMethod);

[data,pos] = collect_box_data(X,box_pos);

centers = mean(box_pos,2);
set(ax,'XTick',centers)
set(ax,'XTickLabel',string(x_list))
xlim(ax,[0.5 nX*nMethod+0.5])

if isempty(data)
    text(ax,0.5,0.5,'No valid estimates', ...
        'Units','normalized','HorizontalAlignment','center');
    draw_group_separators(ax,nX,nMethod)
    return
end

% Missing method/scenario groups must retain their intended empty slots.
% Otherwise boxplot compresses positions and the method colors can shift.
presentPositions = unique(pos,'sorted');
boxplot(ax,data,pos, ...
    'Symbol','', ...
    'Widths',box_width, ...
    'Positions',presentPositions')

hold(ax,'on')

color_boxplot_objects(ax,colors,nMethod)

set(ax,'XTick',centers)
set(ax,'XTickLabel',string(x_list))
xlim(ax,[0.5 nX*nMethod+0.5])

set_adaptive_ylim(ax,ytick_step,isR2)
draw_group_separators(ax,nX,nMethod)

end


%% ============================================================
% Local function 12: Box positions
%% ============================================================
function box_pos = get_box_positions(nX,nMethod)

box_pos = nan(nX,nMethod);

for ix = 1:nX
    start_pos = (ix-1)*nMethod + 1;
    box_pos(ix,:) = start_pos:(start_pos+nMethod-1);
end

end


%% ============================================================
% Local function 13: Collect boxplot data
%% ============================================================
function [data,pos] = collect_box_data(X,box_pos)

[~,nX,nMethod] = size(X);

data = [];
pos  = [];

for ix = 1:nX
    for im = 1:nMethod
        
        x = X(:,ix,im);
        x = x(isfinite(x));
        
        data = [data; x];
        pos  = [pos; box_pos(ix,im)*ones(numel(x),1)];
        
    end
end

end


%% ============================================================
% Local function 14: Color boxplot objects
%% ============================================================
function color_boxplot_objects(ax,colors,nMethod)

h_box = findobj(ax,'Tag','Box');

for i = 1:numel(h_box)
    
    im = method_index_from_position(h_box(i),nMethod);
    
    patch(ax, ...
        get(h_box(i),'XData'), ...
        get(h_box(i),'YData'), ...
        colors(im,:), ...
        'FaceAlpha',0.20, ...
        'EdgeColor',colors(im,:), ...
        'LineWidth',1.5);
    
end

set_colored_lines(ax,'Median',colors,nMethod,2.0)
set_colored_lines(ax,'Whisker',colors,nMethod,1.2)
set_colored_lines(ax,'Upper Whisker',colors,nMethod,1.2)
set_colored_lines(ax,'Lower Whisker',colors,nMethod,1.2)
set_colored_lines(ax,'Upper Adjacent Value',colors,nMethod,1.2)
set_colored_lines(ax,'Lower Adjacent Value',colors,nMethod,1.2)

end


%% ============================================================
% Local function 15: Color line objects
%% ============================================================
function set_colored_lines(ax,tag,colors,nMethod,lineWidth)

h = findobj(ax,'Tag',tag);

for i = 1:numel(h)
    im = method_index_from_position(h(i),nMethod);
    set(h(i),'Color',colors(im,:),'LineWidth',lineWidth)
end

end


function im = method_index_from_position(h,nMethod)
% Read the actual box/line location instead of depending on object order.
x = get(h,'XData');
position = round((min(x)+max(x))/2);
im = mod(position-1,nMethod)+1;
end


%% ============================================================
% Local function 16: Adaptive y-axis
%% ============================================================
function set_adaptive_ylim(ax,ytick_step,isR2)

hU = findobj(ax,'Tag','Upper Adjacent Value');
hL = findobj(ax,'Tag','Lower Adjacent Value');

yu  = collect_ydata(hU);
ylw = collect_ydata(hL);

if isempty(yu) || isempty(ylw)
    return
end

ymin = min(ylw);
ymax = max(yu);

margin = 0.10*(ymax-ymin+eps);

if ymin >= 0
    ylow = 0;
else
    ylow = ymin - margin;
end

yhigh = ymax + margin;

ylow  = floor(ylow/ytick_step)*ytick_step;
yhigh = ceil(yhigh/ytick_step)*ytick_step;

if isR2
    ylow  = max(0,ylow);
    yhigh = min(1,yhigh);
end

ylim(ax,[ylow yhigh])
yticks(ax,ylow:ytick_step:yhigh)

end


%% ============================================================
% Local function 17: Draw group separators
%% ============================================================
function draw_group_separators(ax,nX,nMethod)

yl = ylim(ax);

for ix = 1:nX-1
    sep = ix*nMethod + 0.5;
    plot(ax,[sep sep],yl,'k:','LineWidth',1)
end

ylim(ax,yl)

end


%% ============================================================
% Local function 18: Collect YData
%% ============================================================
function y = collect_ydata(h)

y = [];

for i = 1:numel(h)
    y = [y; h(i).YData(:)];
end

end


%% ============================================================
% Local function 19: Add legend
%% ============================================================
function add_method_legend(ax,method,colors,location)

hold(ax,'on')

h = gobjects(numel(method),1);

for i = 1:numel(method)
    h(i) = plot(ax,nan,nan,'s', ...
        'MarkerFaceColor',colors(i,:), ...
        'MarkerEdgeColor',colors(i,:), ...
        'MarkerSize',8);
end

legend(ax,h,method, ...
    'Location',location, ...
    'Box','on')

end


%% ============================================================
% Local function 20: Axis formatting
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
% Local function 21: Panel label
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
% Local function 22: Save figure
%% ============================================================
function save_figure(fig,filename)

outDir = fullfile(pwd,'figures');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

exportgraphics(fig,fullfile(outDir,filename),'Resolution',600);

end
