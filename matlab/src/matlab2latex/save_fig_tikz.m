function save_fig_tikz(fname,figurehandle,extraopt)
%
% call: save_fig_tikz(fname,figurehandle,extraopt)
%
%
%
if nargin<2 || isempty(figurehandle)
    figurehandle=gcf;
end

if nargin<3 || isempty(extraopt)
    % extraopt='xticklabel style={/pgf/number format/fixed},yticklabel style={/pgf/number format/fixed}';
    extraopt='';
    h=gca;
    objs=h.findobj;
    if isa(objs(end),'matlab.graphics.primitive.Image')
        %extraopt=[extraopt, 'xticklabels={',char(join(string(h.XTickLabel(:,:)),',')),'}'];
        %extraopt=[extraopt, ', yticklabels={',char(join(string(h.YTickLabel(:,:)),',')),'}'];
    else
        extraopt='scaled ticks=false, xticklabel style={/pgf/number format/fixed},yticklabel style={/pgf/number format/fixed}';
    end
end
matlab2tikz(convertStringsToChars(fname+".tex"),'height', '\figureheight', 'width',...
    '\figurewidth','showInfo',false,'figurehandle',figurehandle,'floatFormat','%.6e',...
    'extraAxisOptions',extraopt);
print(figurehandle,[fname],'-dpng', '-r600')
end