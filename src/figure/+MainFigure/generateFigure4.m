function generateFigure4()
% GENERATEFIGURE4_REVISED Generate Figure 4 for the normalization paper
%
% Purpose:
%   Creates a multi-panel figure showing:
%   - Firing rate surfaces for different selectivity/normalization conditions
%   - Relative rate change distributions
%   - Standard deviation analysis across normalization indices
%
% Inputs:
%   None (loads data from files)
%
% Outputs:
%   Saves figure as PDF in results folder
%
% Usage Example:
%   generateFigure4_revised();

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
            'width', 10, 'height', 12);
        
        % Prepare figure
        hFig = figure(1);
        clf;
        utils_plot.matchFigureAspectRatio();
        
        % Create panel layout
        DC1 = utils_plot.divideAxes([0.35, 0.35], [0.25, 0.25, 0.28], ...
            [0.12, 0.07, 0.80, 0.85], 0.04, [0.08, 0.14])';
        DC2 = utils_plot.divideAxes([0.35, 0.35, 0.35], [0.25, 0.25, 0.3], ...
            [0.06, 0.07, 0.93, 0.85], [0.05, 0.12], [0.04, 0.15])';
        
        % Arrange panels
        DC = arrangePanels(DC1, DC2);
        
        % Create axes
        AH = createAxes(DC);
        utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
        % Define color scheme
        cols = copper(6);
        deltaCs = [0.1, 0.2, 0.3, 0.4];
        
        % Plot firing rate surfaces
        plotFiringRateSurfaces(AH(1:4), dataPath, cols, deltaCs);
        
        % Plot relative rate distributions
        plotRateDistributions(AH(5:6), dataPath, cols, deltaCs);
        
        % Plot standard deviation analysis
        plotStandardDeviation(AH(7), dataPath, cols, deltaCs);
        
        % Save figure
        saveFigure(hFig, resultPath, 'Fig4');
        
    catch ME
        fprintf('Error in generateFigure4_revised: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).line);
        end
        rethrow(ME);
    end
end

function DC = arrangePanels(DC1, DC2)
    % Arrange panel positions
    DC = cell(7, 1);
    DC{1} = DC1{1};
    DC{2} = DC1{2};
    DC{3} = DC1{3};
    DC{4} = DC1{4};
    DC{5} = DC2{7};
    DC{6} = DC2{8};
    DC{7} = DC2{9};
end

function AH = createAxes(DC)
    % Create axes for all panels with labels
    
    labels = {'A1', 'A2', 'A3', 'A4', 'B', 'C', 'D'};
    labelPos1 = [-0.08, 0.03];
    labelPos2 = [-0.02, 0.03];
    
    nPanels = length(DC);
    AH = gobjects(nPanels, 1);
    
    for i = 1:nPanels
        AH(i) = axes('Position', DC{i});
        hold(AH(i), 'on');
        
        if i <= 4
            utils_plot.addFigureLabel(labels{i}, labelPos1);
        else
            utils_plot.addFigureLabel(labels{i}, labelPos2);
        end
        
        set(0, 'DefaultAxesTitleFontWeight', 'normal');
    end
end

function plotFiringRateSurfaces(axes_handles, dataPath, cols, deltaCs)
    % Plot firing rate surfaces for different conditions
    
    % Load data
    MTEdata = safeLoadData(fullfile(dataPath, 'Fig4_fr_con1_con2.mat'), 'MTE');
    MTE1data = safeLoadData(fullfile(dataPath, 'Fig1_Pattern3cond.mat'), 'MTE');
    
    % Calculate normalization and selectivity indices
    MTE1 = MTE1data.MTE;
    MTE = MTEdata.MTE;
    
    if size(MTE1, 1) ~= size(MTE, 1)
        error('MTE dimension mismatch between data files');
    end
    
    norm1 = (MTE1(:,1) + MTE1(:,2)) ./ MTE1(:,3);
    sele1 = (MTE1(:,1) - MTE1(:,2)) ./ (MTE1(:,1) + MTE1(:,2));
    
    % Handle NaN and Inf values
    validIdx = isfinite(norm1) & isfinite(sele1);
    norm1 = norm1(validIdx);
    sele1 = sele1(validIdx);
    MTE = MTE(validIdx, :);
    
    markerLineWidth = 1.5;
    
    % Define conditions
    conditions = {
        struct('seleRange', [0, 0.2], 'normRange', [1, 1.4], ...
               'title', 'weak normalization', 'caxis', [0, 60]),
        struct('seleRange', [0, 0.2], 'normRange', [2.6, 3], ...
               'title', 'strong normalization', 'caxis', [0, 15]),
        struct('seleRange', [0.8, 1], 'normRange', [1, 1.4], ...
               'title', '', 'caxis', [0, 40]),
        struct('seleRange', [0.8, 1], 'normRange', [2.6, 3], ...
               'title', '', 'caxis', [0, 30])
    };
    
    % Y-axis labels
    yLabels = {'weak selectivity', '', 'strong selectivity', ''};
    
    for i = 1:4
        axes(axes_handles(i));
        
        % Select data based on conditions
        inds = sele1 > conditions{i}.seleRange(1) & ...
               sele1 < conditions{i}.seleRange(2) & ...
               norm1 > conditions{i}.normRange(1) & ...
               norm1 < conditions{i}.normRange(2);
        
        if sum(inds) == 0
            warning('No data points for condition %d', i);
            continue;
        end
        
        % Reshape and average firing rates
        fr100con1 = reshape(MTE(inds, :), [], 11, 11);
        fr100con1 = squeeze(mean(fr100con1, 1, 'omitnan'));
        
        % Plot surface
        imagesc((0:0.1:1) + 0.05, (0:0.1:1) + 0.05, fr100con1);
        caxis(conditions{i}.caxis);
        colormap('parula');
        
        % Add markers for delta c conditions
        hold on;
        for j = 1:4
            scatter(0.5 + deltaCs(j) + 0.05, 0.5 - deltaCs(j) + 0.05, ...
                30, cols(j,:), 'x', 'LineWidth', markerLineWidth);
        end
        hold off;
        
        % Set properties
        xlim([0, 1.1]);
        ylim([0, 1.1]);
        set(gca, 'YDir', 'normal');
        set(gca, 'DataAspectRatio', [1, 1, 1]);
        
        % Add colorbar
        h = utils_plot.createNarrowColorbar();
        pos = get(h, 'Position');
        set(h, 'Position', pos + [-0.03, 0, 0, 0], 'YAxisLocation', 'right');
        set(h, 'FontSize', 8);
        legend off;
        box off;
        
        % Set tick labels
        set(gca, 'XTick', [0, 0.5, 1] + 0.05);
        set(gca, 'XTickLabels', {'0', '0.5', '1'});
        set(gca, 'YTick', [0, 0.5, 1] + 0.05);
        set(gca, 'YTickLabels', {'0', '0.5', '1'});
        
        % Add labels
        if mod(i, 2) == 1
            ylabel('contrast 2');
            text(-0.5, 0.5, yLabels{i}, 'Units', 'normalized', ...
                'Rotation', 90, 'HorizontalAlignment', 'center', ...
                'FontSize', 10);
        end
        
        if i >= 3
            xlabel('contrast 1');
        end
        
        if i <= 2 && ~isempty(conditions{i}.title)
            text(0.55, 1.3, conditions{i}.title, 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
        
        % Add colorbar label for specific panels
        if i == 2 || i == 4
            text(1.38, 0.5, 'firing rate (Hz)', 'Units', 'normalized', ...
                'Rotation', 90, 'HorizontalAlignment', 'center', ...
                'FontSize', 8);
        end
    end
end

function plotRateDistributions(axes_handles, dataPath, cols, deltaCs)
    % Plot relative rate change distributions
    
    % Load data
    MTEdata = safeLoadData(fullfile(dataPath, 'Fig4_fr_con1_con2.mat'), 'MTE');
    MTE1data = safeLoadData(fullfile(dataPath, 'Fig1_Pattern3cond.mat'), 'MTE');
    
    MTE = MTEdata.MTE;
    MTE1 = MTE1data.MTE;
    
    % Calculate normalized rates
    MTE2 = MTE ./ MTE(:, 61);
    ind1 = find(MTE(:, 61) > 1);
    
    if isempty(ind1)
        warning('No valid data points for rate distributions');
        return;
    end
    
    norm1 = (MTE1(:,1) + MTE1(:,2)) ./ MTE1(:,3);
    norm = norm1(ind1);
    
    % Handle NaN values
    validIdx = isfinite(norm);
    norm = norm(validIdx);
    MTE2 = MTE2(ind1(validIdx), :);
    
    % Setup histogram parameters
    histBins = 0:0.05:3;
    binWidth = 0.05;
    plotX = 0.025:0.05:2.975;
    
    % Condition indices for delta c values
    condInd = [21, 31, 41, 51];
    
    % Plot weak normalization
    axes(axes_handles(1));
    inds = norm >= 1 & norm < 1.4;
    plotDistribution(MTE2(inds, :), condInd, histBins, binWidth, plotX, cols);
    
    pbaspect([1, 1, 1]);
    xlabel('relative rate\newlinechanges (r/r_0)');
    ylabel('probability density');
    
    % Add title
    text('Units', 'Normalized', 'Position', [0.5, 1.08], ...
        'String', 'weak', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
    text('Units', 'Normalized', 'Position', [0.5, 0.95], ...
        'String', 'normalization', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
    
    % Add legend
    addDeltaCLegend(cols, deltaCs);
    
    % Plot strong normalization
    axes(axes_handles(2));
    inds = norm >= 2;
    plotDistribution(MTE2(inds, :), condInd, histBins, binWidth, plotX, cols);
    
    pbaspect([1, 1, 1]);
    ylim([0, 3]);
    xlabel('relative rate\newlinechanges (r/r_0)');
    
    % Add title
    text('Units', 'Normalized', 'Position', [0.5, 1.08], ...
        'String', 'strong', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
    text('Units', 'Normalized', 'Position', [0.5, 0.95], ...
        'String', 'normalization', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
end

function plotStandardDeviation(ax, dataPath, cols, deltaCs)
    % Plot standard deviation across normalization indices
    
    axes(ax);
    
    % Load data
    MTEdata = safeLoadData(fullfile(dataPath, 'Fig4_fr_con1_con2.mat'), 'MTE');
    MTE1data = safeLoadData(fullfile(dataPath, 'Fig1_Pattern3cond.mat'), 'MTE');
    
    MTE = MTEdata.MTE;
    MTE1 = MTE1data.MTE;
    
    % Calculate normalized rates
    MTE2 = MTE ./ MTE(:, 61);
    ind1 = find(MTE(:, 61) > 1);
    
    if isempty(ind1)
        warning('No valid data points for standard deviation analysis');
        return;
    end
    
    norm1 = (MTE1(:,1) + MTE1(:,2)) ./ MTE1(:,3);
    norm = norm1(ind1);
    
    % Handle NaN values
    validIdx = isfinite(norm);
    norm = norm(validIdx);
    MTE2 = MTE2(ind1(validIdx), :);
    
    condInd = [21, 31, 41, 51];
    
    % Calculate standard deviations
    myStd = zeros(length(condInd), 10);
    for ii = 1:10
        inds = norm >= 0.2*ii + 0.8 & norm < 0.2*ii + 1;
        if sum(inds) > 0
            for jj = 1:length(condInd)
                temp1 = MTE2(inds, condInd(jj));
                myStd(jj, ii) = std(temp1, 'omitnan');
            end
        end
    end
    
    % Plot
    hold on;
    for i = 1:4
        plot(1.1:0.2:2.9, myStd(i, :), 'Color', cols(5-i, :), 'LineWidth', 1);
    end
    hold off;
    
    pbaspect([1, 1, 1]);
    xlabel('norm ind');
    ylabel('STD. of r/r0');
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

function plotDistribution(data, condInd, histBins, binWidth, plotX, cols)
    % Plot distribution curves for multiple conditions
    
    densityCurve = zeros(length(condInd), length(plotX));
    
    for jj = 1:length(condInd)
        temp1 = data(:, condInd(jj));
        temp1 = temp1(isfinite(temp1)); % Remove NaN and Inf
        
        if ~isempty(temp1)
            counts = histcounts(temp1, histBins);
            densityCurve(jj, :) = counts / length(temp1) / binWidth;
        end
    end
    
    % Plot curves
    hold on;
    for jj = 1:4
        plot(plotX, densityCurve(jj, :), 'Color', cols(5-jj, :), 'LineWidth', 1);
    end
    hold off;
end

function addDeltaCLegend(cols, deltaCs)
    % Add legend for delta c values
    
    leftLegend = 0.72;
    
    % Add equations
    text('Units', 'Normalized', 'Position', [leftLegend-0.05, 0.8], ...
        'String', '$$c_1 = 0.5+\Delta c$$', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'Color', 'k');
    
    text('Units', 'Normalized', 'Position', [leftLegend-0.05, 0.68], ...
        'String', '$$c_2 = 0.5-\Delta c$$', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'Color', 'k');
    
    text('Units', 'Normalized', 'Position', [leftLegend, 0.56], ...
        'String', '$$\Delta c$$', 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'Color', 'k');
    
    % Add delta c values
    for i = 1:4
        text('Units', 'Normalized', 'Position', [leftLegend, 0.56-0.11*i], ...
            'String', num2str(deltaCs(i)), 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'Color', cols(i, :));
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
