function files=Draw_JoH_Graphical_Abstract(outputDir)
% Reproduce the submitted graphical abstract's layout, palette and visual flow.
% Wording follows the revised manuscript: box fit, then residual PSD feasibility.
% Conceptual curves are illustrative, not additional experimental results.
% Requires base MATLAB only. Exports editable FIG, 300-dpi PNG and vector PDF.
if nargin<1 || isempty(outputDir)
    outputDir=fullfile(fileparts(mfilename('fullpath')),'results','graphical_abstract');
end
if ~isfolder(outputDir), mkdir(outputDir); end
f=figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[2 2 44 23],'DefaultTextFontName','Arial');
guard=onCleanup(@()close(f));
ax=axes(f,'Position',[0 0 1 1]);
axis(ax,[0 100 0 60]); axis(ax,'off'); hold(ax,'on');

blue=[.86 .93 1]; blue2=[.78 .90 1]; red=[1 .87 .87];
yellow=[1 .94 .78]; green=[.88 .97 .88]; grey=[.96 .96 .96];
ecBlue=[.30 .55 .85]; ecRed=[.85 .35 .35];
ecYellow=[.85 .60 .15]; ecGreen=[.35 .70 .40]; ecGrey=[.75 .75 .75];
dark=[.05 .10 .20]; muted=[.35 .35 .35];

label(ax,50,57.6,'N-Tuple Collocation (NC)',25,dark,true);
label(ax,50,55.0,'Robust Reference-Free Error Covariance Estimation',18,dark,true);
h=label(ax,50,52.8,'Beyond Triple Collocation (TC) and Three-Cornered Hat (TCH)',14,muted,false);
h.FontAngle='italic';

% Input and covariance model: same two left-hand panels as the submitted image.
box(ax,4,36,22,14,blue,ecBlue,1.8);
label(ax,15,47.8,'Input',17,dark,true);
label(ax,15,45.2,'Collocated geoscientific datasets',12,dark,true);
label(ax,15,43.0,'No reliable reference truth',12,muted,false);
curves(ax,7,38.3,16,4.1);
label(ax,15,37.2,'x_1, x_2, ..., x_p',14,dark,true);

box(ax,4,20,22,11,grey,ecGrey,1.4);
label(ax,15,28.6,'Covariance-level problem',14,dark,true);
equation(ax,15,25.4,'$\mathbf{Q}=\sigma_y^2\mathbf{a}\mathbf{a}^{\mathrm{T}}+\mathbf{\Sigma}_e$',17,ecBlue);
label(ax,15,22.4,'Estimate \Sigma_e without y',13,muted,false);

% Benchmark panels retain their original colors, positions and short messages.
box(ax,32,38,16,10,red,ecRed,1.8);
label(ax,40,45.8,'TC',19,ecRed,true);
label(ax,40,43.5,'Exactly 3 datasets',12,dark,false);
label(ax,40,41.5,'Uncorrelated errors',12,dark,false);
label(ax,40,39.5,'Negative variance risk',12,ecRed,true);
box(ax,32,24,16,10,yellow,ecYellow,1.8);
label(ax,40,31.8,'TCH',19,ecYellow,true);
label(ax,40,29.5,'Pairwise differences',12,dark,false);
label(ax,40,27.5,'Requires signal cancellation',12,dark,false);
label(ax,40,25.5,'Scaling-sensitive',12,ecYellow,true);

% NC: retain the main blue panel and inset; distinguish fit from feasibility.
box(ax,56,17,24,32,blue2,ecBlue,2.2);
label(ax,68,46,'Proposed NC',20,ecBlue,true);
label(ax,68,43.2,'Constrained covariance decomposition',12,dark,true);
equation(ax,68,40.1,'$\mathbf{Q}=\mathbf{\theta}\mathbf{\theta}^{\mathrm{T}}+\mathbf{\Sigma}_e$',19,ecBlue);
label(ax,68,36.9,'Fit a rank-one signal component',12.5,dark,false);
label(ax,68,34.9,'by minimizing off-diagonal residuals',12.5,dark,false);
equation(ax,68,31.8,'$\min_{\mathbf{\theta}}\ \sum_{i\ne j}|Q_{ij}-\theta_i\theta_j|$',17,ecBlue);
equation(ax,68,28.7,'$\mathrm{subject\ to}\quad\theta_i^2\leq Q_{ii}$',16,ecBlue);
box(ax,57.6,19.2,20.8,6.3,grey,[.70 .82 .95],1.1);
label(ax,68,23.7,'Then enforce residual PSD feasibility',12,dark,true);
label(ax,68,21.2,'Applications: three or more products',12,muted,false);

% Outputs: retain three green panels; avoid claiming universal accuracy.
box(ax,84,39,13,9,green,ecGreen,1.8);
label(ax,90.5,45.9,'Output 1',11.5,muted,false);
label(ax,90.5,43.7,'Error covariance',13,dark,true);
label(ax,90.5,41.5,'without reference',12,dark,false);
box(ax,84,27,13,9,green,ecGreen,1.8);
label(ax,90.5,33.9,'Output 2',11.5,muted,false);
label(ax,90.5,31.7,'Feasible estimates',13,dark,true);
label(ax,90.5,29.5,'PSD residual covariance',11.5,dark,false);
box(ax,84,15,13,9,green,ecGreen,1.8);
label(ax,90.5,21.9,'Output 3',11.5,muted,false);
label(ax,90.5,19.7,'Data merging',13,dark,true);
label(ax,90.5,17.5,'Conditioned WA',12,dark,false);

% Fixed-size arrowheads keep the long lower connector visually restrained.
arrow(f,26,43,31.7,43,dark);
arrow(f,26,25.5,31.7,29,dark);
arrow(f,48.2,43,55.7,40,dark);
arrow(f,48.2,29,55.7,32,dark);
arrow(f,80,43.5,83.7,43.5,dark);
arrow(f,80,31.5,83.7,31.5,dark);
arrow(f,80,19.5,83.7,19.5,dark);
arrow(f,15,20,55.7,20,muted);

box(ax,4,4.5,93,6.7,grey,[.80 .85 .92],1.2);
label(ax,50.5,8.9,'Reference-free uncertainty characterization and multi-source data merging',15.5,dark,true);
label(ax,50.5,6.6,'Assumes signal-uncorrelated errors; physical error covariance is not uniquely identified',13,muted,false);

drawnow;
stem=fullfile(outputDir,'Graphical_Abstract_NC');
files={[stem '.png'],[stem '.pdf'],[stem '.fig']};
exportgraphics(f,files{1},'Resolution',300,'BackgroundColor','white','Padding',12);
exportgraphics(f,files{2},'ContentType','vector','BackgroundColor','white','Padding',12);
savefig(f,files{3});
end

function h=label(ax,x,y,s,sz,c,bold)
h=text(ax,x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Arial','FontSize',sz,'Color',c,'Interpreter','tex');
if bold, h.FontWeight='bold'; end
end

function h=equation(ax,x,y,s,sz,c)
h=text(ax,x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',sz,'Color',c,'Interpreter','latex');
end

function box(ax,x,y,w,h,fc,ec,lw)
rectangle(ax,'Position',[x y w h],'FaceColor',fc,'EdgeColor',ec,'LineWidth',lw);
end

function arrow(f,x1,y1,x2,y2,c)
annotation(f,'arrow',[x1 x2]/100,[y1 y2]/60,'Color',c,'LineWidth',1.6, ...
    'HeadLength',7,'HeadWidth',6);
end

function curves(ax,x,y,w,h)
t=linspace(0,1,120);xx=x+w*t;
plot(ax,xx,y+h*(.72+.12*sin(4*pi*t)),'LineWidth',1.8,'Color',[.80 .10 .10]);
plot(ax,xx,y+h*(.50+.12*sin(4*pi*t+.7)),'LineWidth',1.8,'Color',[0 .35 .75]);
plot(ax,xx,y+h*(.28+.10*sin(4*pi*t-.5)),'LineWidth',1.8,'Color',[0 .55 .20]);
end
