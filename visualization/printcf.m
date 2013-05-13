function printcf(outputfname, fontsize, width, height)
% printcf(outputfname, fontsize, width, height)
%
% Prints the current figure w/ the given dimensions and all text at the
% specified fontsize to a pdf file named in outputfname
%
fontweight = 'Normal';
fontname = 'Times';

% handle subplots appropriately
set(findobj(get(gca, 'Parent'), 'Type', 'axes'), 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname)
titlehs = get(findobj(get(get(gca, 'Parent'), 'Children'), 'Type','axes'), 'Title');
xlabelhs = get(findobj(get(get(gca, 'Parent'), 'Children'), 'Type', 'axes'), 'XLabel');
ylabelhs = get(findobj(get(get(gca, 'Parent'), 'Children'), 'Type', 'axes'), 'YLabel');

% if length(titlehs) > 1
%     for hidx = 1:length(titlehs)
%         set(titlehs{hidx}, 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname);
%     end
% else
%     set(titlehs, 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname);
% end

% if length(xlabelhs) > 1
%     for hidx = 1:length(xlabelhs)
%         set(xlabelhs{hidx}, 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname);
%     end
% else
%     set(xlabelhs, 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname);
% end

 if length(ylabelhs) > 1
     for hidx = 1:length(ylabelhs)
         set(ylabelhs{hidx}, 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname);
     end
 else
     set(ylabelhs, 'FontSize', fontsize, 'FontWeight', fontweight, 'FontName', fontname);
 end

% set(findall(gcf,'type','text'),'FontSize',fontsize)
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperSize', [width height])
set(gcf, 'PaperPosition', [0 0 width height])
print(gcf, '-dpdf', outputfname);
end

