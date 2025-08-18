function generateFigure6()
% GENERATEFIGURE6_REVISED Generate Figure 6 for the normalization paper
%
% Purpose:
%   Creates a multi-panel figure comparing default and in-degree matched conditions:
%   - Normalization index distributions
%   - Fisher Information scaling for contrast and orientation
%   - Neural manifold properties (capacity, radius, dimension)
%
% Inputs:
%   None (loads data from files)
%
% Outputs:
%   Saves figure as PDF in results folder
%
% Usage Example:
%   generateFigure6_revised();

    % Initialize
    try
        % Load constants
        C = figureConstants();
        
        % Setup paths
        dataPath = fullfile(pwd, C.paths.dataFolder);
        resultPath = fullfile(pwd, C.paths.resultFolder);
        
        % Create result directory if needed
        if ~exist(resultPath, 'dir')
            mkdir(resultPath);
        end
        
        % Setup plot options
        utils_plot.setPlotStyle('custom', 'path', resultPath, ...
            'width', 14, 'height', 10);
        
        % Prepare figure
        hFig = figure(1);
        clf;
        utils_plot.matchFigureAspectRatio();
        
        % Create panel layout
        DC1 = utils_plot.divideAxes([0.25, 0.25, 0.25], [0.25, 0.25], ...
            [0.08, 0.09, 0.90, 0.80], [0.12, 0.08], 0.08)';
        
        % Rearrange panels
        DC = cell(6, 1);
        for i = 1:6
            DC{i} = DC1{i};
        end
        
        % Create axes
        AH = createAxes(DC);
        utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
        % Define colors for conditions
        colors = [0, 0, 0; 0.7, 0.7, 0.7];  % Black for default, gray for in-degree
        
        % Plot normalization index distribution
        plotNormalizationDistribution(AH(1), dataPath, colors);
        
        % Plot Fisher Information for contrast
        plotFisherInfoContrast(AH(2), dataPath, colors);
        
        % Plot Fisher Information for orientation
        plotFisherInfoOrientation(AH(3), dataPath, colors);
        
        % Plot neural manifold properties
        plotNeuralManifolds(AH(4:6), dataPath, colors);
        
        % Save figure
        saveFigure(hFig, resultPath, 'Fig6');
        
    catch ME
        fprintf('Error in generateFigure6_revised: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).line);
        end
        rethrow(ME);
    end
end

function AH = createAxes(DC)
    % Create axes for all panels with labels
    
    labels = {'A', 'B', 'C', 'D1', 'D2', 'D3'};
    labelPos = [-0.04, 0.04];
    
    nPanels = length(DC);
    AH = gobjects(nPanels, 1);
    
    for i = 1:nPanels
        AH(i) = axes('Position', DC{i});
        hold(AH(i), 'on');
        utils_plot.addFigureLabel(labels{i}, labelPos);
        set(0, 'DefaultAxesTitleFontWeight', 'normal');
    end
end

function plotNormalizationDistribution(ax, dataPath, colors)
    % Plot normalization index distributions
    
    axes(ax);
    
    % Load data
    defaultData = safeLoadData(fullfile(dataPath, 'Fig2_Pattern3cond.mat'), 'MTE');
    inDegreeData = safeLoadData(fullfile(dataPath, 'FRpattern_InDegree_RightRXScale.mat'), 'MTE');
    
    % Calculate normalization indices
    [normDefault, normInDegree] = calculateNormIndices(defaultData.MTE, inDegreeData.MTE);
    
    % Filter outliers (within 1 std)
    ind11 = filterOutliers(defaultData.MTE);
    ind12 = filterOutliers(inDegreeData.MTE);
    
    normDefault = normDefault(ind11);
    normInDegree = normInDegree(ind12);
    
    % Calculate density curves
    binEdges = 0:0.1:3;
    binCenters = 0.05:0.1:2.95;
    binWidth = 0.1;
    
    densityDefault = calculateDensity(normDefault(normDefault < 3), binEdges, binWidth);
    densityInDegree = calculateDensity(normInDegree(normInDegree < 3), binEdges, binWidth);
    
    % Plot
    hold on;
    plot(binCenters, densityDefault, 'Color', colors(1,:), 'LineWidth', 1);
    plot(binCenters, densityInDegree, 'Color', colors(2,:), 'LineWidth', 1);
    hold off;
    
    xlim([0, 3]);
    ylim([0, 3.5]);
    xlabel('normalization index');
    ylabel('probability density');
    pbaspect([1, 1, 1]);
end

function plotFisherInfoContrast(ax, dataPath, colors)
    % Plot Fisher Information scaling for contrast
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, 'FI_Con.mat'), 'FI_Con');
    
    % Validate data structure
    if ~iscell(data.FI_Con) || length(data.FI_Con) < 2
        error('Invalid FI_Con data structure');
    end
    
    % Setup parameters
    params = struct(...
        'markerSize', 4, ...
        'symbols', {'o', 'square'}, ...
        'nRange', 15000, ...
        'shift', 4000, ...
        'nNeuron', [1,2,4,8,16,31,62,125,250,500,1000,2000,4000,8000,12000], ...
        'upperScale', 1.1, ...
        'nMin', 62);
    
    indPlot = find(params.nNeuron >= params.nMin);
    
    % Plot data
    hold on;
    xlim([params.nMin, params.nRange * 3]);
    
    for i = 1:2
        if isfield(data.FI_Con{i}, 'FI_BC') && isfield(data.FI_Con{i}, 'I')
            % Plot measured points
            plot(params.nNeuron(indPlot), data.FI_Con{i}.FI_BC(indPlot), ...
                params.symbols{i}, 'MarkerSize', params.markerSize, ...
                'Color', colors(i,:));
            
            % Plot asymptotic value
            plot(params.nRange * 3 - params.shift * i, data.FI_Con{i}.I, ...
                params.symbols{i}, 'MarkerSize', params.markerSize, ...
                'Color', colors(i,:));
            
            % Plot theoretical curve
            if isfield(data.FI_Con{i}, 'c')
                theoreticalFI = 1 ./ (1 ./ (data.FI_Con{i}.c * params.nNeuron(indPlot)) + ...
                    1 / data.FI_Con{i}.I);
                plot(params.nNeuron(indPlot), theoreticalFI, ...
                    'Color', colors(i,:), 'LineWidth', 1);
            end
            
            % Add error bars for asymptotic value
            if isfield(data.FI_Con{i}, 'confbound')
                errorbar(params.nRange * 3 - params.shift * i, data.FI_Con{i}.I, ...
                    data.FI_Con{i}.I - data.FI_Con{i}.confbound(1,2), ...
                    data.FI_Con{i}.confbound(2,2) - data.FI_Con{i}.I, ...
                    'Color', colors(i,:));
            end
        end
    end
    
    % Add reference line
    if isfield(data.FI_Con{1}, 'FIin')
        plot(xlim, data.FI_Con{1}.FIin * [1, 1], 'k--', 'LineWidth', 1);
    end
    
    hold off;
    
    % Set properties
    set(gca, 'XScale', 'log');
    ylim([0, data.FI_Con{1}.FIin * params.upperScale]);
    set(gca, 'XTick', [1e2, 1e3, 1e4, params.nRange * 3]);
    set(gca, 'XTickLabel', {'10^2', '10^3', '10^4', '\infty'});
    set(gca, 'YTick', [0, data.FI_Con{1}.FIin/2, data.FI_Con{1}.FIin]);
    set(gca, 'YTickLabel', {'0', 'FI_{in}/2', 'FI_{in}'});
    xlabel('number of neurons');
    ylabel('Fisher information');
    title('Contrast');
    
    % Add legend
    addConditionLegend(colors, {'Default', 'Match in-degree'}, [0.75, 0.35], 0.12);
end

function plotFisherInfoOrientation(ax, dataPath, colors)
    % Plot Fisher Information scaling for orientation
    
    axes(ax);
    
    % Load data
    data = safeLoadData(fullfile(dataPath, 'FI_Ori.mat'), 'FI_Ori');
    
    % Validate data structure
    if ~iscell(data.FI_Ori) || length(data.FI_Ori) < 2
        error('Invalid FI_Ori data structure');
    end
    
    % Setup parameters
    params = struct(...
        'markerSize', 4, ...
        'symbols', {'o', 'square'}, ...
        'nRange', 15000, ...
        'shift', 4000, ...
        'nNeuron', [1,2,4,8,16,31,62,125,250,500,1000,2000,4000,8000,12000], ...
        'upperScale', 1.1, ...
        'nMin', 62);
    
    indPlot = find(params.nNeuron >= params.nMin);
    
    % Plot data
    hold on;
    xlim([params.nMin, params.nRange * 3]);
    
    for i = 1:2
        if isfield(data.FI_Ori{i}, 'FI_BC') && isfield(data.FI_Ori{i}, 'I')
            % Plot measured points
            plot(params.nNeuron(indPlot), data.FI_Ori{i}.FI_BC(indPlot), ...
                params.symbols{i}, 'MarkerSize', params.markerSize, ...
                'Color', colors(i,:));
            
            % Plot asymptotic value
            plot(params.nRange * 3 - params.shift * i, data.FI_Ori{i}.I, ...
                params.symbols{i}, 'MarkerSize', params.markerSize, ...
                'Color', colors(i,:));
            
            % Plot theoretical curve
            if isfield(data.FI_Ori{i}, 'c')
                theoreticalFI = 1 ./ (1 ./ (data.FI_Ori{i}.c * params.nNeuron(indPlot)) + ...
                    1 / data.FI_Ori{i}.I);
                plot(params.nNeuron(indPlot), theoreticalFI, ...
                    'Color', colors(i,:), 'LineWidth', 1);
            end
            
            % Add error bars
            if isfield(data.FI_Ori{i}, 'confbound')
                errorbar(params.nRange * 3 - params.shift * i, data.FI_Ori{i}.I, ...
                    data.FI_Ori{i}.I - data.FI_Ori{i}.confbound(1,2), ...
                    data.FI_Ori{i}.confbound(2,2) - data.FI_Ori{i}.I, ...
                    'Color', colors(i,:));
            end
        end
    end
    
    % Add reference line
    if isfield(data.FI_Ori{1}, 'FIin')
        plot(xlim, data.FI_Ori{1}.FIin * [1, 1], 'k--', 'LineWidth', 1);
    end
    
    hold off;
    
    % Set properties
    set(gca, 'XScale', 'log');
    ylim([0, data.FI_Ori{1}.FIin * params.upperScale]);
    set(gca, 'XTick', [1e2, 1e3, 1e4, params.nRange * 3]);
    set(gca, 'XTickLabel', {'10^2', '10^3', '10^4', '\infty'});
    set(gca, 'YTick', [0, data.FI_Ori{1}.FIin/2, data.FI_Ori{1}.FIin]);
    set(gca, 'YTickLabel', {});
    xlabel('number of neurons');
    title('Orientation');
end

function plotNeuralManifolds(axes_handles, dataPath, colors)
    % Plot neural manifold properties
    
    % Load data
    defaultData = safeLoadData(fullfile(dataPath, 'Neural_Manifolds_Default_50image.mat'), ...
        {'avg_radius_summary', 'avg_dimension_summary', 'avg_capacity_summary'});
    inDegreeData = safeLoadData(fullfile(dataPath, 'Neural_Manifolds_InDegree_50image.mat'), ...
        {'avg_radius_summary', 'avg_dimension_summary', 'avg_capacity_summary'});
    
    % Setup parameters
    nSample = 20;
    nNeuron = [50, 100, 200, 400, 800, 1600, 3200, 6400];
    xLim = [40, 1e4];
    
    % Calculate statistics
    stats = calculateManifoldStats(defaultData, inDegreeData, nSample);
    
    % Plot capacity
    axes(axes_handles(1));
    plotManifoldProperty(stats.capacityMean, stats.capacityStd, ...
        nNeuron, colors, xLim, [0, 0.25], ...
        'capacity', 'number of neurons');
    
    % Plot radius
    axes(axes_handles(2));
    plotManifoldProperty(stats.radiusMean, stats.radiusStd, ...
        nNeuron, colors, xLim, [0.5, 1.1], ...
        'radius', 'number of neurons');
    
    % Plot dimension
    axes(axes_handles(3));
    plotManifoldProperty(stats.dimensionMean, stats.dimensionStd, ...
        nNeuron, colors, xLim, [0, 50], ...
        'dimension', 'number of neurons');
end

% Helper functions
function [normDefault, normInDegree] = calculateNormIndices(MTEdefault, MTEinDegree)
    % Calculate normalization indices
    
    % Validate input dimensions
    if size(MTEdefault, 2) < 3 || size(MTEinDegree, 2) < 3
        error('MTE data must have at least 3 columns');
    end
    
    % Calculate with error handling for division by zero
    normDefault = (MTEdefault(:,1) + MTEdefault(:,2)) ./ MTEdefault(:,3);
    normInDegree = (MTEinDegree(:,1) + MTEinDegree(:,2)) ./ MTEinDegree(:,3);
    
    % Remove invalid values
    normDefault = normDefault(isfinite(normDefault));
    normInDegree = normInDegree(isfinite(normInDegree));
end

function ind = filterOutliers(MTE)
    % Filter data within 1 standard deviation
    
    if size(MTE, 2) < 3
        ind = false(size(MTE, 1), 1);
        return;
    end
    
    ind = true(size(MTE, 1), 1);
    for col = 1:3
        meanVal = mean(MTE(:, col), 'omitnan');
        stdVal = std(MTE(:, col), 'omitnan');
        ind = ind & (abs(MTE(:, col) - meanVal) < stdVal);
    end
end

function density = calculateDensity(data, binEdges, binWidth)
    % Calculate probability density
    
    if isempty(data)
        density = zeros(1, length(binEdges) - 1);
        return;
    end
    
    N = histcounts(data, binEdges);
    density = N / length(data) / binWidth;
end

function stats = calculateManifoldStats(defaultData, inDegreeData, nSample)
    % Calculate manifold statistics
    
    % Initialize
    stats = struct();
    
    % Radius statistics
    stats.radiusMean = [mean(defaultData.avg_radius_summary, 1); ...
                        mean(inDegreeData.avg_radius_summary, 1)];
    stats.radiusStd = [sqrt(var(defaultData.avg_radius_summary, 1) / nSample); ...
                       sqrt(var(inDegreeData.avg_radius_summary, 1) / nSample)];
    
    % Dimension statistics
    stats.dimensionMean = [mean(defaultData.avg_dimension_summary, 1); ...
                          mean(inDegreeData.avg_dimension_summary, 1)];
    stats.dimensionStd = [sqrt(var(defaultData.avg_dimension_summary, 1) / nSample); ...
                         sqrt(var(inDegreeData.avg_dimension_summary, 1) / nSample)];
    
    % Capacity statistics
    stats.capacityMean = [mean(defaultData.avg_capacity_summary, 1); ...
                         mean(inDegreeData.avg_capacity_summary, 1)];
    stats.capacityStd = [sqrt(var(defaultData.avg_capacity_summary, 1) / nSample); ...
                        sqrt(var(inDegreeData.avg_capacity_summary, 1) / nSample)];
end

function plotManifoldProperty(meanData, stdData, nNeuron, colors, xLim, yLim, yLabel, xLabel)
    % Plot manifold property with error bars
    
    % Use only data from 200 neurons onwards (index 3:end)
    plotIdx = 3:length(nNeuron);
    
    for i = 1:2
        if length(meanData(i,:)) >= plotIdx(end) && length(stdData(i,:)) >= plotIdx(end)
            errorbar(nNeuron(plotIdx), meanData(i, plotIdx), ...
                stdData(i, plotIdx), 'Color', colors(i,:), 'LineWidth', 1);
        end
    end
    
    set(gca, 'XScale', 'log');
    xlim(xLim);
    ylim(yLim);
    ylabel(yLabel);
    xlabel(xLabel);
    set(gca, 'XTick', nNeuron);
    pbaspect([1, 1, 1]);
end

function addConditionLegend(colors, labels, startPos, spacing)
    % Add legend for conditions
    
    for i = 1:length(labels)
        text('Units', 'Normalized', ...
            'Position', [startPos(1), startPos(2) - spacing * (i-1)], ...
            'String', labels{i}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 8, ...
            'Color', colors(i,:));
    end
end

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
end.FI_Con{i}.I - data
