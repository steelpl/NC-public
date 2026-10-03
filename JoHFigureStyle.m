function varargout=JoHFigureStyle(action,varargin)
% Shared publication convention. Keep this file identical in both repositories.
% Arial 12 pt axes, 13 pt panel titles, one 11 pt legend, 600 dpi PNG + PDF.
% No global graphics defaults are modified. All drawing/export is MATLAB.
switch lower(action)
    case 'figure'
        w=20; h=12; if ~isempty(varargin), w=varargin{1}; h=varargin{2}; end
        varargout{1}=figure('Visible','off','Color','w','Units','centimeters', ...
            'Position',[2 2 w h],'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
    case 'axes'
        ax=varargin{1}; set(ax,'FontName','Arial','FontSize',12,'LineWidth',.8, ...
            'Box','off','TickDir','out','GridAlpha',.12,'XMinorGrid','off','YMinorGrid','off', ...
            'TitleFontSizeMultiplier',13/12,'LabelFontSizeMultiplier',1);
        grid(ax,'on'); ax.Title.FontWeight='normal';
    case 'legend'
        ax=varargin{1}; tl=varargin{2};
        lg=legend(ax,'show','Interpreter','none','FontName','Arial','FontSize',11, ...
            'Box','off','Orientation','horizontal','NumColumns',3);
        if ~isempty(tl), lg.Layout.Tile='south'; else, lg.Location='southoutside'; end
        varargout{1}=lg;
    case 'color'
        switch string(varargin{1})
            case 'TC', c=[.80 .22 .16];
            case 'TCH', c=[.20 .35 .80];
            case {'NC','NC_PSD'}, c=[.05 .50 .32];
            case 'NC_raw', c=[.80 .55 .05];
            case 'equal_weight', c=[.50 .50 .50];
            otherwise, c=[.12 .12 .12];
        end
        varargout{1}=c;
    case 'export'
        f=varargin{1}; folder=varargin{2}; name=varargin{3};
        if ~isfolder(folder), mkdir(folder); end
        drawnow;
        files={fullfile(folder,[name '.png']),fullfile(folder,[name '.pdf']),fullfile(folder,[name '.fig'])};
        exportgraphics(f,files{1},'Resolution',600,'BackgroundColor','white','Padding',12);
        content='vector'; if numel(varargin)>3, content=varargin{4}; end
        exportgraphics(f,files{2},'ContentType',content,'Resolution',600,'BackgroundColor','white','Padding',12);
        savefig(f,files{3}); varargout{1}=files;
    otherwise
        error('JoHFigure:UnknownAction','Unknown figure-style action: %s',action);
end
end
