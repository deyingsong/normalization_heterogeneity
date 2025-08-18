function SuppFigure_PrefOri_StimOri(varargin)
%SUPPFIGURE_PREFORI_STIMORI Generate supplementary figure for preferred vs stimulus orientation analysis
%
% Purpose:
%   Creates a figure showing correlation patterns between neurons based on
%   their preferred orientation relative to stimulus orientation and absolute
%   preferred orientations, analyzing response correlation structure.
%
% Inputs:
%   varargin - Name-value pairs:
%     'Save'     - logical, whether to save figure (default: true)
%     'View'     - logical, whether to display figure (default: true)
%     'FigNum'   - integer, figure number (default: 1)
%
% Outputs:
%   None (creates and optionally saves figure)
%
% Usage Example:
%   SuppFigure_PrefOri_StimOri();
%   SuppFigure_PrefOri_StimOri('Save', false, 'FigNum', 4);

try
    % Get figure constants
    C = figureConstants();
    thisFile = mfilename('fullpath');
    here     = fileparts(thisFile);          % .../+MainFigure
    rootFolder = fileparts(fileparts(fileparts(here)));  % project root
    resultPath = fullfile(rootFolder, C.paths.resultFolder);
    utilsPath = fullfile(rootFolder, C.paths.utilsFolder);
    figurePath = fullfile(rootFolder, C.paths.figureFolder);
    addpath(resultPath);
    addpath(utilsPath);
    addpath(figurePath);
    
    % Parse input arguments
    p = inputParser;
    addParameter(p, 'Save', true, @islogical);
    addParameter(p, 'View', true, @islogical);
    addParameter(p, 'FigNum', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    
    % Set plot style
    utils_plot.setPlotStyle('custom', 'width', 12, 'height', 5);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = 0.05;
    my_position2 = [0.08 0.18 0.80 0.70];
    
    DC = utils_plot.divideAxes(0.25*ones(1,2), 0.25, my_position2, horizontal_interval, []);
    
    % Load custom colormap
    try
        colormapData = load(fullfile(rootFolder, 'data', 'mycolormap2.mat'));
        if isfield(colormapData, 'mycolormap2')
            customColormap = colormapData.mycolormap2;
        else
            warning('mycolormap2 not found, using default colormap');
            customColormap = parula(256);
        end
    catch ME
        warning('Could not load mycolormap2.mat: %s. Using default colormap.', ME.message);
        customColormap = parula(256);
    end
    
    % Create axes
    LdPos = [-0.03, 0.06];
    labels = {'A', 'B'};
    
    for i = 1:2
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties();
    
    % Parameters
    xticks1 = 0:90:180;
    xticklabels1 = {'-90','0','90'};
    xticklabels2 = {'0', '90', '180'};
    xticks2 = 4.5:9:180;
    
    %% Load data
    try
        rscData = load(fullfile(rootFolder, 'data', 'FigSupp_PrefOri_StimOri_rsc.mat'));
        
        if ~isfield(rscData, 'MTMTpattern')
            error('Required field "MTMTpattern" not found in data file');
        end
        
        MTMTpattern = rscData.MTMTpattern;
        
        % Validate data structure
        if length(MTMTpattern) < 3 || isempty(MTMTpattern{2}) || isempty(MTMTpattern{3})
            error('Insufficient data in MTMTpattern structure');
        end
        
    catch ME
        error('SuppFigure_PrefOri_StimOri:DataLoadError', ...
              'Failed to load required data: %s', ME.message);
    end
    
    % Calculate color limits from all data
    alldata = [MTMTpattern{2}(:); MTMTpattern{3}(:)];
    alldata = alldata(~isnan(alldata) & isfinite(alldata));
    
    if ~isempty(alldata)
        clim1 = round([min(alldata), max(alldata)]/0.01)*0.01;
    else
        warning('No valid data found for color limits, using default');
        clim1 = [-0.05, 0.05];
    end
    
    %% Panel A: Preferred orientation difference vs stimulus orientation
    axes(AH(1));
    
    try
        if ~isempty(MTMTpattern{2}) && ~all(isnan(MTMTpattern{2}(:)))
            imagesc(xticks2, xticks2, MTMTpattern{2}');
            colormap(AH(1), customColormap);
            clim(clim1);
            axis xy;
            
            xlim([0 180]);
            ylim([0 180]);
            xticks(xticks1);
            yticks(xticks1);
            
            xlabel(sprintf('Pref. ori. - stim. ori. (%c)', char(176)), 'FontSize', C.fonts.labelSize);
            ylabel(sprintf('Pref. ori. - stim. ori. (%c)', char(176)), 'FontSize', C.fonts.labelSize);
            xticklabels(xticklabels1);
            yticklabels(xticklabels1);
            
        else
            warning('Panel A data is empty or all NaN');
            imagesc(xticks2, xticks2, rand(length(xticks2), length(xticks2)) * diff(clim1) + clim1(1));
            colormap(AH(1), customColormap);
            clim(clim1);
            axis xy;
            
            text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'Color', 'white', ...
                 'FontSize', C.fonts.annotationSize);
        end
        
    catch ME
        warning('Error processing Panel A: %s', ME.message);
        imagesc(xticks2, xticks2, rand(length(xticks2), length(xticks2)) * diff(clim1) + clim1(1));
        colormap(AH(1), customColormap);
        clim(clim1);
        axis xy;
    end
    
    set(gca, 'DataAspectRatio', [1 1 1]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Absolute preferred orientations
    axes(AH(2));
    
    try
        if ~isempty(MTMTpattern{3}) && ~all(isnan(MTMTpattern{3}(:)))
            imagesc(xticks2, xticks2, MTMTpattern{3}');
            colormap(AH(2), customColormap);
            clim(clim1);
            axis xy;
            
            xlim([0 180]);
            ylim([0 180]);
            xticks(xticks1);
            yticks(xticks1);
            
            xlabel(sprintf('Pref. ori. (%c)', char(176)), 'FontSize', C.fonts.labelSize);
            ylabel(sprintf('Pref. ori. (%c)', char(176)), 'FontSize', C.fonts.labelSize);
            xticklabels(xticklabels2);
            yticklabels(xticklabels2);
            
        else
            warning('Panel B data is empty or all NaN');
            imagesc(xticks2, xticks2, rand(length(xticks2), length(xticks2)) * diff(clim1) + clim1(1));
            colormap(AH(2), customColormap);
            clim(clim1);
            axis xy;
            
            text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'Color', 'white', ...
                 'FontSize', C.fonts.annotationSize);
        end
        
        % Add colorbar for Panel B
        h2 = colorbar;
        set(h2, 'Position', [0.87 0.18 0.013 0.7]);
        clim(clim1);
        
        text('Units', 'Normalized', 'Position', [1.35 0.5], ...
             'string', 'correlation', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, ...
             'Rotation', 90);
        
    catch ME
        warning('Error processing Panel B: %s', ME.message);
        imagesc(xticks2, xticks2, rand(length(xticks2), length(xticks2)) * diff(clim1) + clim1(1));
        colormap(AH(2), customColormap);
        clim(clim1);
        axis xy;
        
        % Still add colorbar even if data failed
        h2 = colorbar;
        set(h2, 'Position', [0.87 0.18 0.013 0.7]);
        clim(clim1);
    end
    
    set(gca, 'DataAspectRatio', [1 1 1]);
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_PrefOri_StimOri', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_PrefOri_StimOri:GeneralError', ...
          'Failed to generate preferred vs stimulus orientation supplementary figure: %s', ME.message);
end

end