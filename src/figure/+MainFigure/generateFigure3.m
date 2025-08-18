function generateFigure3()
% GENERATEFIGURE3_REVISED Generate Figure 3 for the normalization paper
%
% Purpose:
%   Creates a multi-panel figure showing:
%   - Normalization index vs current relationships
%   - Current covariance analyses
%   - Excitatory and inhibitory contributions
%
% Inputs:
%   None (loads data from files)
%
% Outputs:
%   Saves figure as PDF in results folder
%
% Usage Example:
%   generateFigure3_revised();

    % Initialize
    try
        % Load constants
        C = figureConstants();
        
        % Setup paths
        thisFile = mfilename('fullpath');
        here     = fileparts(thisFile);          % .../+MainFigure
        rootFolder = fileparts(fileparts(fileparts(here)));  % project root
        dataPath = fullfile(rootFolder, C.paths.dataFolder);
        resultPath = fullfile(rootFolder, C.paths.resultFolder);
        
        % Create result directory if needed
        if ~exist(resultPath, 'dir')
            mkdir(resultPath);
        end
        
        % Load custom colormap
        colormapData = safeLoadData(fullfile(dataPath, 'mycolormap2.mat'), 'mycolormap2');
        
        % Setup plot options
        utils_plot.setPlotStyle('custom', 'path', resultPath, ...
            'width', 11, 'height', 13.9);
        
        % Prepare figure
        hFig = figure(1);
        clf;
        utils_plot.matchFigureAspectRatio();
        
        % Create panel layout
        DC = utils_plot.divideAxes([0.25, 0.25, 0.25], ...
            [0.25, 0.25, 0.25], [0.07, 0.07, 0.90, 0.85], ...
            0.03 * [1, 1], 0.22 * [1, 1])';
        
        % Create axes
        AH = createAxes(DC);
        utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
        % Define color scheme for Figure 3
        colors = C.colors.fig3_palette;
        
        % Plot correlation scatter plots (top row)
        plotCorrelationScatters(AH(1:3), dataPath, colors, colormapData);
        
        % Plot covariance with fixed axis (middle row)
        plotCovarianceFixedAxis(AH(4:6), dataPath, colors);
        
        % Plot diagonal covariance (bottom row)
        plotCovarianceDiagonal(AH(7:9), dataPath, colors);
        
        % Save figure
        saveFigure(hFig, resultPath, 'Fig3');
        
    catch ME
        fprintf('Error in generateFigure3_revised: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).line);
        end
        rethrow(ME);
    end
end

function AH = createAxes(DC)
    % Create axes for all panels with labels
    
    labels = {'A', 'B', 'C', 'D1', 'D2', 'D3', 'E1', 'E2', 'E3'};
    labelPos = [-0.04, 0.04];
    labelPosBottom = [-0.04, 0.04];
    labelPosMiddle = [0.01, -0.005];
    
    nPanels = numel(DC);
    AH = gobjects(nPanels, 1);
    
    for i = 1:nPanels
        AH(i) = axes('Position', DC{i});
        hold(AH(i), 'on');
        
        if i <= 3
            utils_plot.addFigureLabel(labels{i}, labelPos + [0, 0.005]);
        else
            utils_plot.addFigureLabel(labels{i}, labelPosMiddle);
        end
        
        set(0, 'DefaultAxesTitleFontWeight', 'normal');
    end
end

function plotCorrelationScatters(axes_handles, dataPath, colors, colormapData)
    % Plot normalization index correlation scatter plots
    
    % Load data
    data = safeLoadData(fullfile(dataPath, 'Fig3_NormIndFiringRateCurrent.mat'), ...
        {'normfr', 'normX', 'normE', 'normI', 'rX', 'rE', 'rI'});
    
    % Validate data
    validateVectorSizes({data.normfr, data.normX, data.normE, data.normI});
    
    % Setup common parameters
    NBins = 60;
    xBinWidth = 6 / NBins;
    yBinWidth = 0.3 / NBins;
    sizePoint = 3;
    c0 = parula(256);
    
    % Calculate maximum count for consistent scaling
    N_temp = [histcounts2(data.normfr, data.normX, 0:xBinWidth:6, 1.4:yBinWidth:1.7), ...
              histcounts2(data.normfr, data.normE, 0:xBinWidth:6, 1.4:yBinWidth:1.7), ...
              histcounts2(data.normfr, data.normI, 0:xBinWidth:6, 1.4:yBinWidth:1.7)];
    Nmax = max(N_temp(:));
    
    if Nmax == 0
        warning('No data points in specified range for correlation plots');
        Nmax = 1; % Prevent division by zero
    end
    
    % Data arrays
    dataArrays = {data.normX, data.normE, data.normI};
    correlations = {data.rX, data.rE, data.rI};
    titles = {'Feedforward\newlineexcitation', 'Recurrent\newlineexcitation', 'Recurrent\newlineinhibition'};
    titleColors = {colors(1,:), colors(1,:), colors(5,:)};
    
    for i = 1:3
        axes(axes_handles(i));
        
        % Create 2D histogram
        Ncounts = histcounts2(data.normfr, dataArrays{i}, ...
            0:xBinWidth:6, 1.4:yBinWidth:1.7);
        Ncounts = reshape(Ncounts, [], 1);
        Ninds = ceil(Ncounts / Nmax * 256);
        
        % Create coordinate grid
        [xcoor, ycoor] = meshgrid((0 + xBinWidth/2):xBinWidth:6, ...
            (1.4 + yBinWidth/2):yBinWidth:1.7);
        xcoor = reshape(xcoor', [], 1);
        ycoor = reshape(ycoor', [], 1);
        
        % Create color array
        c = zeros(length(Ninds), 3);
        validInds = Ninds > 0;
        c(validInds, :) = c0(Ninds(validInds), :);
        
        % Remove zero counts
        xcoor = xcoor(validInds);
        ycoor = ycoor(validInds);
        c = c(validInds, :);
        
        % Plot scatter
        if ~isempty(xcoor)
            scatter(xcoor, ycoor, sizePoint, c, 'filled');
        end
        
        % Add correlation value
        if length(correlations{i}) >= 2
            corrValue = correlations{i}(1, 2);
            text('Units', 'Normalized', 'Position', [0.72, 0.89], ...
                'String', sprintf('R=%.2f', corrValue), ...
                'Interpreter', 'latex', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'Color', 'k');
        end
        
        % Set properties
        xlabel('rate norm index');
        if i == 1
            ylabel('current norm index');
        end
        
        % Add title
        text('Units', 'Normalized', 'Position', [0.5, 1.25], ...
            'String', titles{i}(1:end-10), ... % First line
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 9, 'Color', titleColors{i});
        text('Units', 'Normalized', 'Position', [0.5, 1.12], ...
            'String', titles{i}(end-9:end), ... % Second line
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 9, 'Color', titleColors{i});
        
        pbaspect([1, 1, 1]);
        ylim([1.4, 1.7]);
        xlim([0, 6]);
        
        % Add colorbar to last panel
        if i == 3
            h = utils_plot.createNarrowColorbar('vert');
            legend off;
            pos = get(h, 'Position');
            set(h, 'Position', pos + [-0.02, 0, 0, 0]);
            set(h, 'FontSize', 8);
            set(h, 'YTick', [1, 256]);
            set(h, 'YTickLabel', [0, Nmax]);
            
            text(1.15, 0.5, '# neurons', 'Units', 'normalized', ...
                'Rotation', 90, 'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'Color', 'k');
        end
    end
end

function plotCovarianceFixedAxis(axes_handles, dataPath, colors)
    % Plot covariance with one axis fixed
    
    % Load data
    data = safeLoadData(fullfile(dataPath, 'Fig3_CurrentCovarianceFixOneAxis.mat'), ...
        {'covave', 'coverr'});
    
    % Validate data structure
    if ~iscell(data.covave) || ~iscell(data.coverr) || ...
       length(data.covave) < 5 || length(data.coverr) < 5
        error('Invalid covariance data structure');
    end
    
    % Define positions for labels
    positions = struct(...
        'E1E2', [0.8, 0.15], ...
        'I1I2', [0.8, 0.45], ...
        'I1E2', [0.75, 0.7], ...
        'E1I2', [0.75, 0.3], ...
        'total', [0.55, 1.07]);
    
    % Panel 1: E-E and I-I covariance
    axes(axes_handles(1));
    hold on;
    plotErrorbar(1.1:0.2:2.9, data.covave{1}, data.coverr{1}, colors(1,:));
    plotErrorbar(1.1:0.2:2.9, data.covave{4}, data.coverr{4}, colors(5,:));
    hold off;
    ylim([0.1, 0.13]);
    
    addCovarianceLabel('$$\mathrm{Cov}(E_1,E_2)$$', positions.E1E2, colors(1,:));
    addCovarianceLabel('$$\mathrm{Cov}(I_1,I_2)$$', positions.I1I2, colors(5,:));
    
    xlabel('norm index 2');
    ylabel('covariance');
    pbaspect([1, 1, 1]);
    
    % Panel 2: E-I and I-E covariance
    axes(axes_handles(2));
    hold on;
    plotErrorbar(1.1:0.2:2.9, data.covave{2}, data.coverr{2}, colors(4,:));
    plotErrorbar(1.1:0.2:2.9, data.covave{3}, data.coverr{3}, colors(2,:));
    hold off;
    ylim([-0.13, -0.1]);
    
    addCovarianceLabel('$$\mathrm{Cov}(E_1,I_2)$$', positions.E1I2, colors(4,:));
    addCovarianceLabel('$$\mathrm{Cov}(I_1,E_2)$$', positions.I1E2, colors(2,:));
    
    xlabel('norm index 2');
    title('1 < norm index 1 < 1.5', 'Units', 'normalized', 'Position', [0.5, 1.3]);
    pbaspect([1, 1, 1]);
    
    % Panel 3: Total covariance
    axes(axes_handles(3));
    plotErrorbar(1.1:0.2:2.9, data.covave{5}, data.coverr{5}, 'k');
    
    addCovarianceLabel('$$\mathrm{Cov}(\mathrm{Total}_1,\mathrm{Total}_2)$$', positions.total, 'k');
    
    xlabel('norm index 2');
    pbaspect([1, 1, 1]);
end

function plotCovarianceDiagonal(axes_handles, dataPath, colors)
    % Plot diagonal covariance analysis
    
    % Load data
    data = safeLoadData(fullfile(dataPath, 'Fig3_CurrentCovarianceDiagonal.mat'), ...
        {'covave', 'coverr'});
    
    % Validate data structure
    if ~iscell(data.covave) || ~iscell(data.coverr) || ...
       length(data.covave) < 5 || length(data.coverr) < 5
        error('Invalid diagonal covariance data structure');
    end
    
    % Define positions for labels
    positions = struct(...
        'E1E2', [0.8, 0.15], ...
        'I1I2', [0.8, 0.45], ...
        'I1E2_2', [0.4, 0.58], ...
        'E1I2_2', [0.4, 0.39], ...
        'total', [0.55, 1.07]);
    
    % Panel 1: E-E and I-I covariance
    axes(axes_handles(1));
    hold on;
    plotErrorbar(1.1:0.2:2.9, data.covave{1}, data.coverr{1}, colors(1,:));
    plotErrorbar(1.1:0.2:2.9, data.covave{4}, data.coverr{4}, colors(5,:));
    hold off;
    ylim([0.1, 0.14]);
    
    addCovarianceLabel('$$\mathrm{Cov}(E_1,E_2)$$', positions.E1E2, colors(1,:));
    addCovarianceLabel('$$\mathrm{Cov}(I_1,I_2)$$', positions.I1I2, colors(5,:));
    
    xlabel('norm index');
    ylabel('covariance');
    pbaspect([1, 1, 1]);
    
    % Panel 2: E-I covariance
    axes(axes_handles(2));
    plotErrorbar(1.1:0.2:2.9, data.covave{2}, data.coverr{2}, colors(4,:));
    ylim([-0.14, -0.1]);
    
    addCovarianceLabel('$$\mathrm{Cov}(E_1,I_2)$$', positions.I1E2_2, colors(4,:));
    addCovarianceLabel('$$(=\mathrm{Cov}(E_2,I_1))$$', positions.E1I2_2, colors(4,:));
    
    xlabel('norm index');
    title('similar norm index', 'Units', 'normalized', 'Position', [0.5, 1.3]);
    pbaspect([1, 1, 1]);
    
    % Panel 3: Total covariance
    axes(axes_handles(3));
    plotErrorbar(1.1:0.2:2.9, data.covave{5}, data.coverr{5}, 'k');
    
    addCovarianceLabel('$$\mathrm{Cov}(\mathrm{Total}_1,\mathrm{Total}_2)$$', positions.total, 'k');
    
    xlabel('norm index');
    pbaspect([1, 1, 1]);
end

% Helper functions
function data = safeLoadData(filepath, requiredFields)
    % Load data with error handling and field validation
    
    if ~exist(filepath, 'file')
        error('Data file not found: %s', filepath);
    end
    
    try
        data = load(filepath);
    catch ME
        error('Failed to load %s: %s', filepath, ME.message);
    end
    
    % Validate required fields
    if nargin > 1
        if ischar(requiredFields)
            requiredFields = {requiredFields};
        end
        
        for i = 1:length(requiredFields)
            if ~isfield(data, requiredFields{i})
                error('Required field "%s" not found in %s', ...
                    requiredFields{i}, filepath);
            end
        end
    end
end

function validateVectorSizes(vectors)
    % Validate that all vectors have the same size
    
    if isempty(vectors)
        return;
    end
    
    refSize = size(vectors{1});
    for i = 2:length(vectors)
        if ~isequal(size(vectors{i}), refSize)
            error('Vector size mismatch: expected %s, got %s', ...
                mat2str(refSize), mat2str(size(vectors{i})));
        end
    end
end

function plotErrorbar(x, y, err, color)
    % Plot errorbar with consistent styling
    
    if isempty(y) || isempty(err)
        warning('Empty data for errorbar plot');
        return;
    end
    
    % Handle NaN values
    validIdx = ~isnan(y) & ~isnan(err);
    if sum(validIdx) < 2
        warning('Insufficient valid data points for errorbar plot');
        return;
    end
    
    errorbar(x(validIdx), y(validIdx), err(validIdx), ...
        'Color', color, 'LineWidth', 1);
end

function addCovarianceLabel(labelText, position, color)
    % Add covariance label with LaTeX formatting
    
    text('Units', 'Normalized', 'Position', position, ...
        'String', labelText, ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8, ...
        'Color', color);
end

function saveFigure(hFig, resultPath, filename)
    % Save figure with error handling
    
    try
        set(hFig, 'Renderer', 'painters');
        utils_plot.saveFigure('path', resultPath, 'name', filename, ...
            'view', 1, 'save', 1, 'format', 'pdf', 'res', 600);
        fprintf('Figure saved successfully: %s.pdf\n', filename);
    catch ME
        warning('Failed to save figure: %s', ME.message);
    end
end