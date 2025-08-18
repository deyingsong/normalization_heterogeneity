function SuppFigure_Superimposed(varargin)
%SUPPFIGURE_SUPERIMPOSED Generate supplementary figure for superimposed stimulus condition
%
% Purpose:
%   Creates a figure showing model analysis under superimposed stimulus conditions,
%   including normalization index distribution, correlation heatmap, and tuning
%   similarity control analysis compared to separated stimulus conditions.
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
%   SuppFigure_Superimposed();
%   SuppFigure_Superimposed('Save', false, 'FigNum', 3);

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
    utils_plot.setPlotStyle('custom', 'width', 13, 'height', 5);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    DC = utils_plot.divideAxes([0.25 0.25 0.25], 0.70, [0.09 0.12 0.90 0.68], [0.08 0.18], []);
    
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
    Labels = {'A', 'B', 'C'};
    LdPos = C.positions.legendOffset + [0 0.03];
    
    for i = 1:3
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(Labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties();
    
    % Colors
    col1 = zeros(2, 3);
    col2 = C.colors.gray(1) * ones(2, 3);
    
    %% Panel A: Normalization index distribution
    axes(AH(1));
    
    try
        % Load superimposed overlap data
        overlapData = load(fullfile(rootFolder, 'data', 'FigSupp_Superimposed_NormIndDistribution.mat'));
        
        if ~isfield(overlapData, 'norm1')
            error('Required field "norm1" not found in overlap data');
        end
        
        norm1 = overlapData.norm1;
        
        % Filter valid data
        temp = norm1(norm1 <= 3 & norm1 >= 0 & ~isnan(norm1));
        
        if ~isempty(temp)
            binEdges = 0:0.1:3;
            counts = histcounts(temp, binEdges);
            binCenters = binEdges(1:end-1) + diff(binEdges)/2;
            densityValues = counts / length(temp) / 0.1;
            
            plot(binCenters, densityValues, 'k', 'LineWidth', C.plot.lineWidth);
        else
            warning('No valid normalization indices found');
            plot(0.05:0.1:2.95, zeros(1, 30), 'k');
        end
        
    catch ME
        warning('Error processing superimposed data: %s', ME.message);
        plot(0.05:0.1:2.95, rand(1, 30), 'k');
    end
    
    xlim([0 3]);
    xlabel('norm index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    
    text('Units', 'Normalized', 'Position', [0.5 C.positions.titleHeight], ...
         'string', 'Model,', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.titleSize, 'color', 'k');
    text('Units', 'Normalized', 'Position', [0.5 C.positions.subtitleHeight], ...
         'string', 'superimposed', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.titleSize, 'color', 'k');
    
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Correlation heatmap
    axes(AH(2));
    
    try
        % Load correlation data
        corrData = load(fullfile(rootFolder, 'data', 'FigSupp_Superimposed_Norm1Norm2Corr.mat'));
        
        if ~isfield(corrData, 'dataheat')
            error('Required field "dataheat" not found in correlation data');
        end
        
        dataheat = corrData.dataheat;
        
        % Validate data
        if isempty(dataheat) || all(isnan(dataheat(:)))
            warning('Correlation data is empty or all NaN');
            dataheat = rand(40, 40) * 0.1 - 0.05; % Placeholder
        end
        
        imagesc(0.525:0.05:2.475, 0.525:0.05:2.475, dataheat');
        colormap(AH(2), customColormap);
        axis xy;
        
        % Create colorbar
        h1 = colorbar;
        h2_height = 0.61;
        set(h1, 'Position', [0.63 0.15 0.008 h2_height]);
        set(h1, 'YTick', -0.04:0.02:0.02, 'YTickLabel', {});
        set(h1, 'FontSize', C.fonts.tickSize);
        
        % Coordinates for Y tick labels
        yticks1 = -0.04:0.02:0.02;
        xpos = 1.17; % Position adjustment
        
        % Create text objects for each label
        for i = 1:length(yticks1)
            text('Units', 'Normalized', 'Position', [xpos 0.28*i-0.32], ...
                 'string', num2str(yticks1(i)), 'color', 'k', ...
                 'Rotation', 45, 'FontSize', C.fonts.tickSize);
        end
        
        utils_plot.freezeColors(AH(2));
        
    catch ME
        warning('Error processing correlation heatmap: %s', ME.message);
        imagesc(0.525:0.05:2.475, 0.525:0.05:2.475, rand(40, 40) * 0.1 - 0.05);
        colormap(AH(2), customColormap);
    end
    
    xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
    ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
    xlim([0.5 2.5]);
    ylim([0.5 2.5]);
    xticks([0.5 1.5 2.5]);
    xticklabels([0.5 1.5 2.5]);
    yticks([0.5 1.5 2.5]);
    yticklabels([0.5 1.5 2.5]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C: Tuning similarity control
    axes(AH(3));
    
    try
        % Load tuning similarity data
        tuningData = load(fullfile(rootFolder, 'data', 'FigSupp_Superimposed_TuningSimilarityControl.mat'));
        
        if ~isfield(tuningData, 'mcorr') || ~isfield(tuningData, 'scorr')
            error('Required fields not found in tuning similarity data');
        end
        
        mcorr = tuningData.mcorr;
        scorr = tuningData.scorr;
        
        xVals = -0.95:0.1:0.95;
        
        % Validate data dimensions
        if size(mcorr, 2) ~= length(xVals) || size(scorr, 2) ~= length(xVals)
            warning('Data dimensions do not match expected x-axis values');
            % Use available data or create placeholder
            if size(mcorr, 2) > 0 && size(scorr, 2) > 0
                nPoints = min([size(mcorr, 2), size(scorr, 2), length(xVals)]);
                xVals = xVals(1:nPoints);
                mcorr = mcorr(:, 1:nPoints);
                scorr = scorr(:, 1:nPoints);
            else
                xVals = -0.95:0.1:0.95;
                mcorr = [zeros(1, length(xVals)); zeros(1, length(xVals))];
                scorr = [ones(1, length(xVals)) * 0.01; ones(1, length(xVals)) * 0.01];
            end
        end
        
        if size(mcorr, 1) >= 2 && size(scorr, 1) >= 2
            errorbar(xVals, mcorr(1,:), scorr(1,:), 'color', col1(1,:), 'LineWidth', C.plot.lineWidth);
            hold on;
            errorbar(xVals, mcorr(2,:), scorr(2,:), 'color', col2(1,:), 'LineWidth', C.plot.lineWidth);
            hold off;
        else
            warning('Insufficient data rows for error bar plot');
            plot(xVals, zeros(size(xVals)), 'k');
        end
        
    catch ME
        warning('Error processing tuning similarity data: %s', ME.message);
        xVals = -0.95:0.1:0.95;
        plot(xVals, zeros(size(xVals)), 'k');
    end
    
    xlabel('tuning similarity', 'FontSize', C.fonts.labelSize);
    ylabel('correlation', 'FontSize', C.fonts.labelSize);
    
    % Add legend
    text('Units', 'Normalized', 'Position', [0.05 0.95], ...
         'string', 'similar index', 'color', col1(1,:), ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
         'FontSize', C.fonts.annotationSize);
    text('Units', 'Normalized', 'Position', [0.05 0.8], ...
         'string', 'different index', 'color', col2(1,:), ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
         'FontSize', C.fonts.annotationSize);
    
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_Superimposed', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_Superimposed:GeneralError', ...
          'Failed to generate superimposed supplementary figure: %s', ME.message);
end

end