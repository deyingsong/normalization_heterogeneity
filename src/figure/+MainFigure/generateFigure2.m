function generateFigure2()
% GENERATEFIGURE2_REVISED Generate Figure 2 for the normalization paper
%
% Purpose:
%   Creates a multi-panel figure comparing model and experimental data:
%   - Normalization index correlation heat maps
%   - Diagonal and anti-diagonal correlation analyses  
%   - Tuning similarity controls
%
% Inputs:
%   None (loads data from files)
%
% Outputs:
%   Saves figure as PDF in results folder
%
% Usage Example:
%   generateFigure2_revised();

    % Initialize
    try
        % Load constants
        C = figureConstants();
        
        % Setup paths
        parentFolder = fileparts(pwd);
        dataPath = fullfile(parentFolder, C.paths.dataFolder);
        resultPath = fullfile(parentFolder, C.paths.resultFolder);
        
        % Create result directory if needed
        if ~exist(resultPath, 'dir')
            mkdir(resultPath);
        end
        
        % Load custom colormap
        colormapData = safeLoadData(fullfile(dataPath, 'mycolormap2.mat'), 'mycolormap2');
        
        % Setup plot options
        utils_plot.setPlotStyle('custom', 'path', resultPath, ...
            'width', 16, 'height', 12);
        
        % Prepare figure
        hFig = figure(1);
        clf;
        utils_plot.matchFigureAspectRatio();
        
        % Create panel layout
        DC = utils_plot.divideAxes([0.22, 0.22, 0.22, 0.22], ...
            [0.4, 0.4, 0.4], [0.06, 0.08, 0.90, 0.86], ...
            [0.10, 0.02, 0.02], [0.2, 0.2])';
        
        % Create axes
        AH = createAxes(DC);
        utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
        % Define colors
        grayColor = C.colors.gray;
        
        % Plot model data (top row)
        plotModelHeatmap(AH(1), dataPath, colormapData);
        plotModelAntiDiagonal(AH(2), dataPath);
        plotModelDiagonal(AH(3), dataPath);
        plotModelTuningSimilarity(AH(4), dataPath, grayColor);
        
        % Plot MT data (middle row)
        plotMTHeatmap(AH(5), dataPath, colormapData);
        plotMTAntiDiagonal(AH(6), dataPath);
        plotMTDiagonal(AH(7), dataPath);
        plotMTTuningSimilarity(AH(8), dataPath, grayColor);
        
        % Plot V4 data (bottom row)
        plotV4Heatmap(AH(9), dataPath, colormapData);
        plotV4AntiDiagonal(AH(10), dataPath);
        plotV4Diagonal(AH(11), dataPath);
        plotV4TuningSimilarity(AH(12), dataPath, grayColor);
        
        % Save figure
        saveFigure(hFig, resultPath, 'Fig2');
        
    catch ME
        fprintf('Error in generateFigure2_revised: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).line);
        end
        rethrow(ME);
    end
end

function AH = createAxes(DC)
    % Create axes for all panels with labels
    
    labels = {'A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'C1', 'C2', 'C3', 'C4'};
    labelPos = [-0.03, 0.03];
    
    nPanels = numel(DC);
    AH = gobjects(nPanels, 1);
    
    for i = 1:nPanels
        AH(i) = axes('Position', DC{i});
        hold(AH(i), 'on');
        utils_plot.addFigureLabel(labels{i}, labelPos);
        set(0, 'DefaultAxesTitleFontWeight', 'normal');
    end
end

function plotModelHeatmap(ax, dataPath, colormapData)
    % Plot model normalization index correlation heatmap
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, 'Fig2_Norm1Norm2Corr.mat'), 'dataheat');
    
    % Plot heatmap
    imagesc(1.025:0.05:2.975, 1.025:0.05:2.975, data.dataheat');
    caxis([-0.01, 0.02]);
    box on;
    
    % Set colormap
    colormap(ax, colormapData.mycolormap2);
    
    % Add colorbar
    h1 = utils_plot.createNarrowColorbar('vert');
    set(h1, 'YLim', [-0.01, 0.02]);
    set(h1, 'YTick', -0.01:0.01:0.02, 'YTickLabel', {}, ...
        'YTickLabelRotation', 45, 'FontSize', 8);
    
    % Add colorbar tick labels manually
    addColorbarLabels(h1, -0.01:0.01:0.02, 7);
    
    % Set properties
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    axis xy;
    xlabel('norm index 1');
    ylabel('norm index 2');
    xlim([1, 3]);
    ylim([1, 3]);
    xticks([1, 2, 3]);
    yticks([1, 2, 3]);
    
    title('Model');
    
    % Add correlation label
    addRotatedLabel('correlation', [1.45, 0.5], 90, 8);
end

function plotModelAntiDiagonal(ax, dataPath)
    % Plot model anti-diagonal correlation analysis
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig2_AntiDiagonalCorrDecreaseWithDiffNMI.mat'), {'ave', 'err'});
    
    % Plot with error bars
    errorbar(0.05:0.1:0.45, data.ave, data.err, 'k', 'LineWidth', 1);
    
    % Set properties
    pbaspect([1, 1, 1]);
    xlabel('index 1 - index 2');
    xlim([0, 0.5]);
end

function plotModelDiagonal(ax, dataPath)
    % Plot model diagonal correlation analysis
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig2_DiagonalCorrDecreaseWithAveNMI.mat'), {'ave', 'err'});
    
    % Plot with error bars
    errorbar(1.05:0.1:2.95, data.ave, data.err, 'k', 'LineWidth', 1);
    
    % Set properties
    pbaspect([1, 1, 1]);
    xlabel('norm index');
    xlim([1, 3]);
end

function plotModelTuningSimilarity(ax, dataPath, grayColor)
    % Plot model tuning similarity control
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig2_TuningSimilarityControl.mat'), {'mcorr', 'scorr'});
    
    % Validate data dimensions
    if size(data.mcorr, 2) < 19 || size(data.scorr, 2) < 19
        error('Invalid tuning similarity data dimensions');
    end
    
    % Plot both conditions
    errorbar(-0.95:0.1:0.95, data.mcorr(1,:), data.scorr(1,:), ...
        'k', 'LineWidth', 1);
    hold on;
    errorbar(-0.95:0.1:0.95, data.mcorr(2,:), data.scorr(2,:), ...
        'Color', grayColor, 'LineWidth', 1);
    hold off;
    
    % Set properties
    xlabel('tuning similarity');
    pbaspect([1, 1, 1]);
    xlim([-1, 1]);
    
    % Add legend text
    text('Units', 'Normalized', 'Position', [0.05, 0.95], ...
        'String', 'similar index', 'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'FontSize', 8);
    text('Units', 'Normalized', 'Position', [0.05, 0.8], ...
        'String', 'different index', 'Color', grayColor, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'FontSize', 8);
end

function plotMTHeatmap(ax, dataPath, colormapData)
    % Plot MT normalization index correlation heatmap
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_MT_Norm1Norm2Corr.mat'), 'dataheat1');
    
    % Plot heatmap
    imagesc(1.05:0.1:1.95, 1.05:0.1:1.95, data.dataheat1');
    box on;
    
    % Set colormap
    colormap(ax, colormapData.mycolormap2);
    
    % Add colorbar
    h1 = utils_plot.createNarrowColorbar('vert');
    set(h1, 'YTick', 0.05:0.05:0.25, 'YTickLabel', {}, ...
        'YTickLabelRotation', 45, 'FontSize', 8);
    
    % Add colorbar tick labels
    addColorbarLabels(h1, 0.05:0.05:0.25, 7);
    
    % Set properties
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    axis xy;
    xlabel('norm index 1');
    ylabel('norm index 2');
    xlim([1, 2]);
    ylim([1, 2]);
    xticks([1, 1.5, 2]);
    yticks([1, 1.5, 2]);
    
    title('MT data');
    
    % Add correlation label
    addRotatedLabel('correlation', [1.5, 0.5], 90, 8);
end

function plotMTAntiDiagonal(ax, dataPath)
    % Plot MT anti-diagonal correlation analysis
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_MT_AntiDiagonalCorrDecreaseWithDiffNMI.mat'), ...
        {'myave', 'myerr'});
    
    % Plot with error bars
    errorbar(0.05:0.1:0.45, data.myave, data.myerr, 'k', 'LineWidth', 1);
    
    % Set properties
    pbaspect([1, 1, 1]);
    xlabel('index 1 - index 2');
    xlim([0, 0.5]);
    yticks(0:0.1:0.3);
end

function plotMTDiagonal(ax, dataPath)
    % Plot MT diagonal correlation analysis
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_MT_DiagonalCorrDecreaseWithAveNMI_Width_0d2.mat'), ...
        {'myave', 'myerr'});
    
    % Plot with error bars
    errorbar(1.05:0.1:1.95, data.myave, data.myerr, 'k', 'LineWidth', 1);
    
    % Set properties
    xlabel('norm index');
    pbaspect([1, 1, 1]);
    xlim([1, 2]);
end

function plotMTTuningSimilarity(ax, dataPath, grayColor)
    % Plot MT tuning similarity control
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_MT_TuningSimilarityControl_weighted.mat'), ...
        {'myaves', 'myerrs', 'myaved', 'myerrd'});
    
    % Plot both conditions
    errorbar(-0.9:0.2:0.9, data.myaves, data.myerrs, ...
        'k', 'LineWidth', 1);
    hold on;
    errorbar(-0.9:0.2:0.9, data.myaved, data.myerrd, ...
        'Color', grayColor, 'LineWidth', 1);
    hold off;
    
    % Set properties
    xlabel('tuning similarity');
    pbaspect([1, 1, 1]);
    xlim([-1, 1]);
    
    % Add legend text
    addLegendText({'similar index', 'different index'}, ...
        {[0.05, 0.95], [0.05, 0.8]}, {'k', grayColor}, 8);
end

function plotV4Heatmap(ax, dataPath, colormapData)
    % Plot V4 normalization index correlation heatmap
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_V4_Norm1Norm2Corr.mat'), 'dataheat1');
    
    % Plot heatmap
    imagesc(1.05:0.1:1.95, 1.05:0.1:1.95, data.dataheat1');
    box on;
    
    % Set colormap
    colormap(ax, colormapData.mycolormap2);
    
    % Add colorbar
    h1 = utils_plot.createNarrowColorbar('vert');
    set(h1, 'YTick', 0.04:0.02:0.12, 'YTickLabel', {}, ...
        'YTickLabelRotation', 45, 'FontSize', 8);
    
    % Add colorbar tick labels
    addColorbarLabels(h1, 0.04:0.02:0.12, 7);
    
    % Set properties
    set(gca, 'DataAspectRatio', [1, 1, 1]);
    axis xy;
    xlabel('norm index 1');
    ylabel('norm index 2');
    xlim([1, 2]);
    ylim([1, 2]);
    xticks([1, 1.5, 2]);
    yticks([1, 1.5, 2]);
    
    title('V4 data');
    
    % Add correlation label
    addRotatedLabel('correlation', [1.5, 0.5], 90, 8);
end

function plotV4AntiDiagonal(ax, dataPath)
    % Plot V4 anti-diagonal correlation analysis
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_V4_AntiDiagonalCorrDecreaseWithDiffNMI.mat'), ...
        {'myave', 'myerr'});
    
    % Plot with error bars
    errorbar(0.05:0.1:0.45, data.myave, data.myerr, 'k', 'LineWidth', 1);
    
    % Set properties
    pbaspect([1, 1, 1]);
    xlabel('index 1 - index 2');
    xlim([0, 0.5]);
end

function plotV4Diagonal(ax, dataPath)
    % Plot V4 diagonal correlation analysis
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_V4_DiagonalCorrDecreaseWithAveNMI_Width_0d2.mat'), ...
        {'myave', 'myerr'});
    
    % Plot with error bars
    errorbar(1.05:0.1:1.95, data.myave, data.myerr, 'k', 'LineWidth', 1);
    
    % Set properties
    xlabel('norm index');
    pbaspect([1, 1, 1]);
    xlim([1, 2]);
end

function plotV4TuningSimilarity(ax, dataPath, grayColor)
    % Plot V4 tuning similarity control
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig_ExpData_V4_TuningSimilarityControl_weighted.mat'), ...
        {'myaves', 'myerrs', 'myaved', 'myerrd'});
    
    % Plot both conditions
    errorbar(-0.9:0.2:0.9, data.myaves, data.myerrs, ...
        'k', 'LineWidth', 1);
    hold on;
    errorbar(-0.9:0.2:0.9, data.myaved, data.myerrd, ...
        'Color', grayColor, 'LineWidth', 1);
    hold off;
    
    % Set properties
    xlabel('tuning similarity');
    pbaspect([1, 1, 1]);
    xlim([-1, 1]);
    
    % Add legend text
    addLegendText({'similar index', 'different index'}, ...
        {[0.05, 0.95], [0.05, 0.8]}, {'k', grayColor}, 8);
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

function addColorbarLabels(h, ticks, xpos)
    % Add manual tick labels to colorbar
    
    for i = 1:length(ticks)
        text(xpos, ticks(i), num2str(ticks(i)), ...
            'Parent', h, 'Rotation', 45, 'FontSize', 8);
    end
end

function addRotatedLabel(labelText, position, rotation, fontSize)
    % Add rotated text label
    
    h = text('Units', 'Normalized', 'Position', position, ...
        'String', labelText, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', fontSize);
    set(h, 'Rotation', rotation);
end

function addLegendText(labels, positions, colors, fontSize)
    % Add legend text annotations
    
    for i = 1:length(labels)
        text('Units', 'Normalized', 'Position', positions{i}, ...
            'String', labels{i}, 'Color', colors{i}, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', fontSize);
    end
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
