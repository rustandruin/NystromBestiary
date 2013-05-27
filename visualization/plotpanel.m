function plotpanel(varargin)
%
% plotpanel(series)
%
% Plots the information in the structure series to the current axis with a
% legend
%
% plotpanel(series, legendflag)
%
% Does the same. Legend is plotted if legendflag is true, otherwise it is
% not
%
% series.plotname -- the name of the plot (optional)
% series.seriesnames -- a string cell array with the legend names for the series
% series.xlabel -- label for the xaxis (optional)
% series.ylabel -- label for the yaxis (optional)
% series.clipaxis -- a flag; true indicates to not let the yaxis have negative values (optional)
% series.legendloc -- location for the legend
% series.plottype -- type of plot (optional)
%
% series.x -- the common x-values each series should be plotted against
% series.series -- a cell array containing the data series (as vectors)
% series.styles -- a string array cell containing the plot styles for each series
% series.lw -- a cell containing the linewidths for the series
% series.ms -- a cell containing the marker sizes for the series
% series.colors -- a cell array containing the rgb colors for the series (as vectors)
% series.mcolors -- a cell array containing the rgb colors for the series (as vectors)
%

series = varargin{1};

x = series.x;
data = series.series;
styles = series.styles;
lw = series.lw;
ms = series.ms;
colors = series.colors;
mcolors = series.mcolors;

if isfield(series, 'plottype')
    switch lower(series.plottype)
        case 'semilogy'
            myplot = @semilogy;
        otherwise
            myplot = @plot;
    end
else
    myplot = @plot;
end

for idx = 1:length(data)
    myplot(x, data{idx}, styles{idx}, 'LineWidth', lw{idx}, ...
        'MarkerSize', ms{idx}, 'Color', colors{idx}, ...
        'MarkerFaceColor', mcolors{idx});
    hold on;
end
hold off;

if (isfield(series, 'clipaxis') && series.clipaxis)
    extents = axis();
    axis([extents(1) extents(2) max(extents(3), 0), extents(4)])
end

if isfield(series, 'xlabel')
    xlabel(series.xlabel, 'Interpreter', 'LaTex', 'fontsize', 20);
end
if isfield(series, 'ylabel')
    ylabel(series.ylabel, 'fontsize', 20);
end

showlegendq = true;
if nargin == 2
 showlegendq = varargin{2};
end
if showlegendq
    legend(char(series.seriesnames), 'Location', series.legendloc);
end

if isfield(series, 'plotname')
    th = title(series.plotname, 'fontsize', 20);
    set(th, 'Interpreter', 'LaTex');
end

end
