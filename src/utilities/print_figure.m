function print_figure(fig_handle,varargin)
% PRINT_FIGURE   Figure printing utility
% This function is useful to generate figures that are uniformly looking in
% terms of fonts and size, for example for a thesis. It can produce
% high-resolution bitmaps (png) or vector graphics (eps with high-res
% preview, or pdf).
%
% PRINT_FIGURE(fig_handle) generates an image file in the current
% directory; the file name is the same as the figure name. All other
% settings are at their default values.
% 
% Optional arguments:
% 
% PRINT_FIGURE(fig_handle, size) resizes the figure to match the size
% specified:
% 's': small (width 3", height 2")
% 'c': column-wide for two-column papers (width 8.5 cm = 3.35 in, height 6 cm = 2.35 in)
% 'm': medium (width 4.5", height 3") 
% 'l': large (width 5.5", height 3")
% 'w': wide (width 6", height 3")
% 'q': square (width 4", height 4")
% [width height]: custom (give width and height in inches)
% 'screen': maintain screen dimension (no resize) [default]
%
% PRINT_FIGURE(fig_handle, size, filename) uses the string specified in filename
% as name for the file. Default filename is "fig_#", where # is the figure
% number.
%
% PRINT_FIGURE(fig_handle, size, filename, format) uses the format specified by
% the string "format": at the moment, this can be:
% 'png' (300 dpi) [default]
% 'pdf'
% 'eps' (with TIFF preview)
% Use png for hi-res bitmap graphic, good for slides; eps or pdf for 
% vectorial graphic that will print neatly, especially in LaTeX. 
% (The eps driver used is epsc2, i.e. color level 2 PostScript, with tiff preview)
%
% PRINT_FIGURE(fig_handle, size, filename, format, fontname, fontsize) forces the
% font face and size of the axis labels, legend, title, and annotations to
% the specified values. Enter a string with the font name and a numerical
% value for the font size.
% 
% All inputs except for the the first are optional, but their order is fixed:
% enter an empty array [] to leave an option at its default value.
%
%
% Default settings:
%    size = 'screen'
%    filename = "Fig_#" where # is the figure handle (number)
%    format = 'png' (300 dpi)
%    fontname = 'Helvetica'
%    fontsize = 10
% 
% Example:
%  figure, plot(10*rand(1,40)), xlabel('time [s]'), ylabel('V [m/s]')
%  print_figure(gcf,'m',[],'png')
%  print_figure(gcf,'c','Driving_cycle','pdf','Times New Roman',8)
% 
%
% See also PRINT

if nargin >= 2 && ~isempty(varargin{1})
    FigureSize = varargin{1};
else
    FigureSize = 'screen'; % default
end

if nargin >= 3 && ~isempty(varargin{2})
    filename = varargin{2};
else
    filename = ['Fig_' num2str(fig_handle)]; % default
end

if nargin >= 4 && ~isempty(varargin{3})
    FileFormat = varargin{3};
else
    FileFormat = 'png'; % default
end

if ~strcmpi(FileFormat,{'png','eps','pdf'})
    error('Supported file formats are PNG, EPS, PDF')
end

if nargin >= 5 && ~isempty(varargin{4})
    font_name = varargin{4};
else
    font_name = 'Helvetica'; % default
end

if nargin >= 6 && ~isempty(varargin{5})
    font_size= varargin{5};
else
    font_size = 10;
end

% Text formatting
hAll = findall(gcf); % get all objects in the figure
for idx = 1 : length(hAll)
    try % try to set font properties, if possible
        set(hAll(idx),'fontname',font_name,'FontSize',font_size);
%         set(hAll(idx),'interpreter','LaTeX')
    catch
        % ... otherwise never mind...
    end
end

%% printing figure
% resize

% check if figure is docked (resize doesn't work for docked figures:
origWindowStyle = get(fig_handle,'WindowStyle');
if strcmpi(origWindowStyle,'docked')
    set(fig_handle,'WindowStyle','normal')
end

% decide new size (using the "position" property)
origUnits = get(fig_handle,'units');
set(fig_handle,'units','inches');
origPosition = get(fig_handle,'position');
newPosition = origPosition;

if isnumeric(FigureSize) & numel(FigureSize) == 2 % "custom" option
    newPosition(3:4) = FigureSize;
elseif ischar(FigureSize) % one of the predefined options
    switch FigureSize
        case 's' % small, 3x2"
            newPosition(3:4) = [3 2];
        case 'c' % column, 3.35x2.35"
%             newPosition(3:4) = [3.35 2.35];
            newPosition(3:4) = [3.35 2.25];
        case 'm' % medium, 4.5x3" [default]
            newPosition(3:4) = [4.5 3];
        case 'l' % large, 5.5x3"
            newPosition(3:4) = [5.5 3];
        case 'w' % wide, 6x3"
            newPosition(3:4) = [6 3];
        case 'q' % square, 4x4"
            newPosition(3:4) = [4 4];
        case 'screen' % no resize, screen dimension
            newPosition(3:4) = origPosition(3:4);
        otherwise
            newPosition(3:4) = origPosition(3:4); % [default]
            warning('The option for size was not recognized. Setting size to screen value.')
    end
else % not a valid input
    error('Invalid option for figure size')
end
set(fig_handle,'position',newPosition)

% print...
set(fig_handle,...
    'PaperUnits','inches',...
    'PaperSize', [newPosition(3) newPosition(4)+0.3],... % give a bit of extra vertical space
    'paperPositionMode','auto')

switch FileFormat
    case 'png'
        print(fig_handle,'-dpng','-r300',filename) % PNG
    case 'eps'
        print(fig_handle, '-depsc2', '-r300', filename) % color EPS with TIFF preview
    case 'pdf'        
        print(fig_handle, '-dpdf', filename) % PDF
        % to embed fonts:
%         eval(['!/usr/local/bin/ps2pdf13 -dPDFSETTINGS=/prepress -dPDFCrop ' filename '.pdf' filename '.pdf'])
%         print_pdf(filename, fig_handle)
end

% ...and then return things as they were
set(fig_handle,'position',origPosition);

set(fig_handle,'units',origUnits)

if strcmpi(origWindowStyle,'docked')
    set(fig_handle,'WindowStyle','docked') 
end

% to embed fonts in pdf, use this system command (pdf conversion utility):
% eval(['!ps2pdf13 -dPDFSETTINGS=/prepress ' filename '.pdf' filename '.pdf'])

