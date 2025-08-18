function matchFigureAspectRatio(varargin)
% MATCHFIGUREASPECTRATIO Match figure window to paper aspect ratio
%
% Syntax:
%   matchFigureAspectRatio()
%   matchFigureAspectRatio('Name', Value, ...)
%
% Description:
%   Adjusts figure window size to match paper dimensions for accurate preview.
%
% Inputs (Name-Value pairs):
%   'Area'   - 'print' (default) or 'keep' - sizing method
%   'Factor' - Scaling factor (default: 1)
%
% Example:
%   matchFigureAspectRatio('Factor', 0.5);

    % Constants
    PIXELS_PER_INCH = get(0, 'ScreenPixelsPerInch');
    CM_TO_INCH = 2.54;
    DEFAULT_PIXELS_PER_CM = 44.9; % Default for typical display
    
    % Parse parameters
    params = utils_plot.parseNameValuePairs(varargin{:});
    params = utils_plot.validateField(params, 'Area', 'print', 'char');
    params = utils_plot.validateField(params, 'Factor', 1, 'numeric');
    
    % Get current figure properties
    paperPos = get(gcf, 'PaperPosition');
    aspectRatio = paperPos(4) / paperPos(3);
    figPos = get(gcf, 'Position');
    
    % Calculate conversion factor
    paperUnits = get(gcf, 'PaperUnits');
    conversionFactor = calculateConversionFactor(paperUnits, PIXELS_PER_INCH, CM_TO_INCH);
    
    % Calculate new figure size
    switch lower(params.Area)
        case 'keep'
            % Maintain figure area, adjust aspect ratio
            figArea = figPos(3) * figPos(4);
            newWidth = round(sqrt(figArea / aspectRatio));
            newHeight = round(newWidth * aspectRatio);
            
        case 'print'
            % Match paper size
            newWidth = round(conversionFactor * paperPos(3) / params.Factor);
            newHeight = round(conversionFactor * paperPos(4) / params.Factor);
            
        otherwise
            error('matchFigureAspectRatio:InvalidArea', ...
                  'Area must be "print" or "keep"');
    end
    
    % Apply new size
    figPos(3:4) = [newWidth, newHeight];
    set(gcf, 'Position', figPos);
    
    % --- Helper Functions ---
    function factor = calculateConversionFactor(units, ppi, cm2in)
        switch lower(units)
            case 'inches'
                factor = ppi;
            case 'centimeters'
                factor = ppi / cm2in;
            otherwise
                warning('matchFigureAspectRatio:UnknownUnits', ...
                        'Unknown paper units, using default conversion');
                factor = DEFAULT_PIXELS_PER_CM;
        end
    end
end
