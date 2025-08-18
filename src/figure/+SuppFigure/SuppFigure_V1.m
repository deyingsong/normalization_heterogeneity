function SuppFigure_V1(varargin)
%SUPPFIGURE_V1 Generate supplementary figure showing V1 experimental data analysis
%
% Purpose:
%   Creates a three-panel figure showing V1 experimental data including
%   normalization index distribution, correlation analysis, and tuning
%   similarity control analysis.
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
%   SuppFigure_V1();                    % Default behavior
%   SuppFigure_V1('Save', false);       % Don't save
%   SuppFigure_V1('FigNum', 5);         % Use figure 5

try
    % Get figure constants
    C = figureConstants();
    thisFile = mfilename('fullpath');
    here     = fileparts(thisFile);          % .../+MainFigure
    rootFolder = fileparts(fileparts(fileparts(here)));  % project root
    resultPath = fullfile(rootFolder, C.paths.resultFolder);
    utilsPath = fullfile(rootFolder, C.paths.utilsFolder);
    addpath(resultPath);
    addpath(utilsPath);
    
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
    DC = utils_plot.divideAxes([0.25 0.25 0.25], 0.70, [0.09 0.12 0.90 0.68], ...
                               [0.08 0.18], []);
    
    % Load custom colormap
    try
        colormapData = load(fullfile(rootFolder,'data', 'mycolormap2.mat'));
        if isfield(colormapData, 'mycolormap2')
            customColormap = colormapData.mycolormap2;
        else
            warning('mycolormap2 not found in file, using default colormap');
            customColormap = parula(256);
        end
    catch ME
        warning('Could not load mycolormap2.mat: %s. Using default colormap.', ME.message);
        customColormap = parula(256);
    end
    
    % Create axes
    Labels = {'A', 'B', 'C'};
    LdPos = C.positions.legendOffset+[0 0.05];
    
    for i = 1:3
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(Labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Colors
    col1 = zeros(2, 3);
    col2 = C.colors.gray(1) * ones(2, 3);
    
    %% Panel A: Normalization index distribution
    axes(AH(1));
    
    try
        % Load V1 normalization index data
        v1Data = load(fullfile(rootFolder,'data', 'Fig_ExpData_V1_NormIndDistribution.mat'));
        
        if ~isfield(v1Data, 'norminds')
            error('Required field "norminds" not found in V1 data file');
        end
        
        norminds = v1Data.norminds;
        
        % Validate data
        if isempty(norminds)
            error('V1 normalization indices data is empty');
        end
        
        % Remove NaN and extreme values
        norminds = norminds(~isnan(norminds) & norminds <= 3 & norminds >= 0);
        
        if isempty(norminds)
            error('No valid normalization indices after filtering');
        end
        
        % Calculate histogram
        binEdges = 0:0.1:3;
        counts = histcounts(norminds, binEdges);
        binCenters = binEdges(1:end-1) + diff(binEdges)/2;
        densityValues = counts / length(norminds) / 0.1;
        
        plot(binCenters, densityValues, 'k', 'LineWidth', C.plot.lineWidth);
        
    catch ME
        warning('Error processing V1 data for Panel A: %s', ME.message);
        % Create placeholder plot
        plot(0.05:0.1:2.95, zeros(1, 30), 'k');
        text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center');
    end
    
    xlim([0 3]);
    xlabel('norm index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    
    text('Units', 'Normalized', 'Position', [0.5 C.positions.subtitleHeight], ...
         'string', 'V1 data', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.titleSize, ...
         'color', 'k');
    
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Correlation heatmap
    axes(AH(2));
    
    try
        % Load correlation data
        corrData = load(fullfile(rootFolder, 'data', 'Fig_ExpData_V1_Norm1Norm2Corr.mat'));
        
        if ~isfield(corrData, 'dataheat1')
            error('Required field "dataheat1" not found in correlation data file');
        end
        
        dataheat = corrData.dataheat1;
        
        % Validate data
        if isempty(dataheat) || any(isnan(dataheat(:)))
            warning('Correlation data contains NaN or empty values');
        end
        
        imagesc(1.05:0.1:1.95, 1.05:0.1:1.95, dataheat');
        colormap(AH(2), customColormap);
        
        % Create colorbar
        h1 = colorbar;
        h2_height = 0.62;
        set(h1, 'Position', [0.62 0.15 0.008 h2_height]);
        set(h1, 'YTick', 0:0.05:0.15, 'YTickLabel', {});
        set(h1, 'FontSize', C.fonts.tickSize);
        
        % Add tick labels manually
        yticks1 = 0:0.05:0.15;
        xpos = 1.17;
        
        for i = 1:length(yticks1)
            text('Units', 'Normalized', 'Position', [xpos 0.31*i-0.33], ...
                 'string', num2str(yticks1(i)), 'color', 'k', ...
                 'Rotation', 45, 'FontSize', C.fonts.tickSize);
        end
        
    catch ME
        warning('Error processing correlation data for Panel B: %s', ME.message);
        % Create placeholder
        imagesc(1.05:0.1:1.95, 1.05:0.1:1.95, rand(9, 9));
        text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'Color', 'white');
    end
    
    axis xy;
    xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
    ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
    xlim([1 2]);
    ylim([1 2]);
    xticks([1 1.5 2]);
    xticklabels([1 1.5 2]);
    yticks([1 1.5 2]);
    yticklabels([1 1.5 2]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C: Tuning similarity control
    axes(AH(3));
    
    try
        % Load tuning similarity data
        tuningData = load(fullfile(rootFolder, 'data', 'Fig_ExpData_V1_TuningSimilarityControl.mat'));
        
        requiredFields = {'myaves', 'myerrs', 'myaved', 'myerrd'};
        for field = requiredFields
            if ~isfield(tuningData, field{1})
                error('Required field "%s" not found in tuning data file', field{1});
            end
        end
        
        xVals = -0.9:0.2:0.9;
        
        % Validate data dimensions
        if length(tuningData.myaves) ~= length(xVals) || ...
           length(tuningData.myerrs) ~= length(xVals) || ...
           length(tuningData.myaved) ~= length(xVals) || ...
           length(tuningData.myerrd) ~= length(xVals)
            error('Data dimensions do not match expected x-axis values');
        end
        
        % Plot error bars
        errorbar(xVals, tuningData.myaves, tuningData.myerrs, ...
                'color', col1(2,:), 'LineWidth', C.plot.lineWidth);
        hold on;
        errorbar(xVals, tuningData.myaved, tuningData.myerrd, ...
                'color', col2(2,:), 'LineWidth', C.plot.lineWidth);
        hold off;
        
    catch ME
        warning('Error processing tuning similarity data for Panel C: %s', ME.message);
        % Create placeholder
        xVals = -0.9:0.2:0.9;
        plot(xVals, zeros(size(xVals)), 'k');
        text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center');
    end
    
    xlabel('tuning similarity', 'FontSize', C.fonts.labelSize);
    ylabel('correlation', 'FontSize', C.fonts.labelSize);
    
    % Add legend
    text('Units', 'Normalized', 'Position', [0.05 0.95], ...
         'string', 'similar index', 'color', col1(2,:), ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
         'FontSize', C.fonts.annotationSize);
    text('Units', 'Normalized', 'Position', [0.05 0.8], ...
         'string', 'different index', 'color', col2(2,:), ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
         'FontSize', C.fonts.annotationSize);
    
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_V1', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        % Just show the figure
        drawnow;
    end
    
catch ME
    error('SuppFigure_V1:GeneralError', ...
          'Failed to generate supplementary figure V1: %s', ME.message);
end

end