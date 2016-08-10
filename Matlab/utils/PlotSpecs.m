function [lineSpecs,textSpecs,plotSpecs] = PlotSpecs()

lineWidth = 2;
markerSize = 10;
fontSize = 12;
fontWeight = 'bold';
fontColor = 'k';

textSpecs = {'FontSize',fontSize,'FontWeight',fontWeight,'Color',fontColor};
lineSpecs = {'LineWidth',lineWidth,'MarkerSize',markerSize};
plotSpecs = {'numbertitle','off','WindowStyle','docked'};

%% Old stuff
% if ~exist('style','var')
%     style = 'presentation';
% end
% 
% if strcmpi(style,'paper')
%     % Defaults for this papers
%     width = 3;     % Width in inches
%     height = 3;    % Height in inches
%     alw = 0.75;    % AxesLineWidth
%     fsz = 8;      % Fontsize
%     lw = 1.5;      % LineWidth
%     msz = 8;       % MarkerSize
% elseif strcmpi(style,'presentation')
%     % Defaults for slides
%     width = 6;     % Width in inches
%     height = 5;    % Height in inches
%     alw = 1;    % AxesLineWidth
%     fsz = 14;      % Fontsize
%     lw = 2;      % LineWidth
%     msz = 12;       % MarkerSize
% else
%     
% end
% 
% 
% % The new defaults will not take effect if there are any open figures. To
% % use them, we close all figures, and then repeat the first example.
% close all;
% 
% % The properties we've been using in the figures
% set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
% set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% 
% % Set the default Size for display
% defpos = get(0,'defaultFigurePosition');
% set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);
% 
% % Set the defaults for saving/printing to a file
% set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
% set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
% defsize = get(gcf, 'PaperSize');
% left = (defsize(1)- width)/2;
% bottom = (defsize(2)- height)/2;
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);
% 
% close all;
end