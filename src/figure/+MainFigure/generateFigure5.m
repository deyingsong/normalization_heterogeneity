function generateFigure5()
% GENERATEFIGURE5_REVISED Generate Figure 5 for the normalization paper
%
% Purpose:
%   Creates a multi-panel figure showing Fisher Information analysis:
%   - FI per spike for contrast discrimination
%   - FI per spike for orientation discrimination
%   - Grouped by normalization index ranges
%
% Inputs:
%   None (loads data from files)
%
% Outputs:
%   Saves figure as PDF in results folder
%
% Usage Example:
%   generateFigure5_revised();

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
            'width', 12, 'height', 9);
        
        % Prepare figure
        hFig = figure(1);
        clf;
        utils_plot.matchFigureAspectRatio();
        
        % Create panel layout
        DC1 = utils_plot.divideAxes([0.25, 0.25, 0.25], [0.25, 0.25], ...
            [0.11, 0.09, 0.85, 0.80], 0.05 * [1, 1], 0.15)';
        
        % Rearrange panels
        DC = cell(6, 1);
        DC{1} = DC1{1};
        DC{2} = DC1{2};
        DC{3} = DC1{3};
        DC{4} = DC1{4};
        DC{5} = DC1{5};
        DC{6} = DC1{6};
        
        % Create axes
        AH = createAxes(DC);
        utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
        % Define parameters
        colors = hot(8);
        idPlots = [1, 5];  % Which normalization groups to highlight
        
        % Process and plot contrast FI data
        plotContrastFI(AH(1:3), dataPath, colors, idPlots);
        
        % Process and plot orientation FI data
        plotOrientationFI(AH(4:6), dataPath, colors, idPlots);
        
        % Save figure
        saveFigure(hFig, resultPath, 'Fig5');
        
    catch ME
        fprintf('Error in generateFigure5_revised: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).line);
        end
        rethrow(ME);
    end
end

function AH = createAxes(DC)
    % Create axes for all panels with labels
    
    labels = {'A1', 'A2', 'A3', 'B1', 'B2', 'B3'};
    labelPos1 = [-0.04, 0.04];
    labelPos2 = [-0.04, 0.04];
    
    nPanels = length(DC);
    AH = gobjects(nPanels, 1);
    
    for i = 1:nPanels
        AH(i) = axes('Position', DC{i});
        hold(AH(i), 'on');
        
        if i <= 3
            utils_plot.addFigureLabel(labels{i}, labelPos2);
        else
            utils_plot.addFigureLabel(labels{i}, labelPos1);
        end
        
        set(0, 'DefaultAxesTitleFontWeight', 'normal');
    end
end

function plotContrastFI(axes_handles, dataPath, colors, idPlots)
    % Plot Fisher Information for contrast discrimination
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig5_FI_Con_Default_RandomNum_Arithmetic.mat'), ...
        'FI_Con_Default_RandomNum_Arithmetic');
    
    % Process data into groups
    [conSpkCnt, conFI] = processFIData(data.FI_Con_Default_RandomNum_Arithmetic);
    
    % Calculate binned statistics
    plotFICon = calculateBinnedStats(conSpkCnt, conFI);
    
    % Define plot parameters
    params = struct(...
        'nRange', 1.1e4, ...
        'nBins', 10, ...
        'markerAlpha', 0.2, ...
        'xTicks', {[1, 1e2, 1e4]}, ...
        'xTickLabels', {{'10^0', '10^2', '10^4'}}, ...
        'yLim', {[0, 10]}, ...
        'yLimAll', {[0, 5]}, ...
        'titleHeight', 1.25, ...
        'labels', {{'[1,1.4)', '[1.4,1.8)', '[1.8,2.2)', '[2.2,2.6)', '[2.6,3)'}});
    
    % Plot highlighted groups
    plotHighlightedGroup(axes_handles(1), conSpkCnt{idPlots(1)}, ...
        conFI{idPlots(1)}, plotFICon{idPlots(1)}, ...
        colors(idPlots(1), :), params, 1, params.labels{1});
    
    plotHighlightedGroup(axes_handles(2), conSpkCnt{idPlots(2)}, ...
        conFI{idPlots(2)}, plotFICon{idPlots(2)}, ...
        colors(idPlots(2), :), params, 2, params.labels{5});
    
    % Add title to middle panel
    axes(axes_handles(2));
    text('Units', 'Normalized', 'Position', [0.5, params.titleHeight], ...
        'String', 'Contrast', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
    
    % Plot all groups
    axes(axes_handles(3));
    plotAllGroups(plotFICon, colors, params);
    
    % Add legend
    addNormIndexLegend(params.labels, colors);
end

function plotOrientationFI(axes_handles, dataPath, colors, idPlots)
    % Plot Fisher Information for orientation discrimination
    
    % Load data
    data = safeLoadData(fullfile(dataPath, ...
        'Fig5_FI_Ori_Default_RandomNum_Arithmetic.mat'), ...
        'FI_Ori_Default_RandomNum_Arithmetic');
    
    % Process data into groups
    [oriSpkCnt, oriFI] = processFIData(data.FI_Ori_Default_RandomNum_Arithmetic);
    
    % Calculate binned statistics
    plotFIOri = calculateBinnedStats(oriSpkCnt, oriFI);
    
    % Define plot parameters
    params = struct(...
        'nRange', 1.1e4, ...
        'nBins', 10, ...
        'markerAlpha', 0.2, ...
        'xTicks', {[1, 1e2, 1e4]}, ...
        'xTickLabels', {{'10^0', '10^2', '10^4'}}, ...
        'yLim', {[0, 1.7]}, ...
        'yLimAll', {[0, 1]}, ...
        'titleHeight', 1.25, ...
        'labels', {{'[1,1.4)', '[1.4,1.8)', '[1.8,2.2)', '[2.2,2.6)', '[2.6,3)'}});
    
    % Plot highlighted groups
    plotHighlightedGroupOri(axes_handles(1), oriSpkCnt{idPlots(1)}, ...
        oriFI{idPlots(1)}, plotFIOri{idPlots(1)}, ...
        colors(idPlots(1), :), params, true);
    
    plotHighlightedGroupOri(axes_handles(2), oriSpkCnt{idPlots(2)}, ...
        oriFI{idPlots(2)}, plotFIOri{idPlots(2)}, ...
        colors(idPlots(2), :), params, false);
    
    % Add title to middle panel
    axes(axes_handles(2));
    text('Units', 'Normalized', 'Position', [0.5, params.titleHeight], ...
        'String', 'Orientation', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k');
    
    % Plot all groups
    axes(axes_handles(3));
    plotAllGroupsOri(plotFIOri, colors, params);
end

% Helper functions
function [spkCnt, FI] = processFIData(dataCell)
    % Process Fisher Information data from cell array structure
    
    nGroups = 5;
    spkCnt = cell(nGroups, 1);
    FI = cell(nGroups, 1);
    
    for i = 1:nGroups
        flag = 0;
        tempSpk = [];
        tempFI = [];
        
        for j = 1:size(dataCell, 2)
            if ~isempty(dataCell{i, j}) && isfield(dataCell{i, j}, 'spkcnt')
                nSamples = dataCell{i, j}.Nsample;
                if nSamples > 0
                    tempSpk = [tempSpk, dataCell{i, j}.spkcnt];
                    tempFI = [tempFI, dataCell{i, j}.FI];
                end
            end
        end
        
        spkCnt{i} = tempSpk;
        FI{i} = tempFI;
    end
end

function plotFI = calculateBinnedStats(spkCnt, FI)
    % Calculate binned statistics for FI data
    
    nGroups = length(spkCnt);
    nBins = 10;
    plotFI = cell(nGroups, 1);
    
    for idNormGroup = 1:nGroups
        if isempty(spkCnt{idNormGroup})
            plotFI{idNormGroup} = zeros(3, 0);
            continue;
        end
        
        spkMax = max(spkCnt{idNormGroup});
        if spkMax <= 0
            plotFI{idNormGroup} = zeros(3, 0);
            continue;
        end
        
        logSpkMax = log(spkMax);
        plotFI{idNormGroup} = zeros(3, nBins);
        
        for i = 1:nBins
            binMin = logSpkMax / nBins * (i - 1);
            binMax = logSpkMax / nBins * i;
            
            inBin = log(spkCnt{idNormGroup}) >= binMin & ...
                    log(spkCnt{idNormGroup}) < binMax;
            
            if sum(inBin) > 0
                tempFI = FI{idNormGroup}(inBin);
                tempSpk = spkCnt{idNormGroup}(inBin);
                normalizedFI = tempFI ./ tempSpk;
                
                % Remove outliers
                validIdx = isfinite(normalizedFI);
                if sum(validIdx) > 0
                    plotFI{idNormGroup}(1, i) = logSpkMax / nBins * (i - 0.5);
                    plotFI{idNormGroup}(2, i) = mean(normalizedFI(validIdx));
                    plotFI{idNormGroup}(3, i) = std(normalizedFI(validIdx)) / ...
                        sqrt(sum(validIdx));
                end
            end
        end
    end
end

function plotHighlightedGroup(ax, spkCnt, FI, plotFI, color, params, panelNum, label)
    % Plot highlighted group for contrast with scatter and errorbar
    axes(ax);
    hold on;
    
    % Scatter plot
    if ~isempty(spkCnt) && ~isempty(FI)
        normalizedFI = FI ./ spkCnt;
        validIdx = isfinite(normalizedFI) & spkCnt > 0;
        
        s = scatter(spkCnt(validIdx), normalizedFI(validIdx), 1, color, 'filled');
        s.MarkerFaceAlpha = params.markerAlpha;
        s.MarkerEdgeAlpha = params.markerAlpha;
    end
    
    % Error bar plot
    if size(plotFI, 2) > 0
        validCols = plotFI(2, :) > 0;
        errorbar(exp(plotFI(1, validCols)), plotFI(2, validCols), ...
            plotFI(3, validCols), 'Color', color);
    end
    
    hold off;
    
    % Set properties
    set(gca, 'XScale', 'log');
    ylim(params.yLim);
    xlim([1, params.nRange]);
    set(gca, 'XTick', params.xTicks);
    set(gca, 'XTickLabel', {});
    
    % Add labels for panel 1
    if panelNum == 1
        addAxisLabels();
        
        % Add norm index label
        text('Units', 'Normalized', 'Position', [0.8, 0.9], ...
            'String', 'norm. ind.', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', 'k');
        text('Units', 'Normalized', 'Position', [0.8, 0.8], ...
            'String', label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', color);
    elseif panelNum == 2
        % Add norm index label
        text('Units', 'Normalized', 'Position', [0.8, 0.9], ...
            'String', 'norm. ind.', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', 'k');
        text('Units', 'Normalized', 'Position', [0.8, 0.8], ...
            'String', label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', color);
    end
end

function plotHighlightedGroupOri(ax, spkCnt, FI, plotFI, color, params, showLabels)
    % Plot highlighted group for orientation
    axes(ax);
    hold on;
    
    % Scatter plot
    if ~isempty(spkCnt) && ~isempty(FI)
        normalizedFI = FI ./ spkCnt;
        validIdx = isfinite(normalizedFI) & spkCnt > 0;
        
        s = scatter(spkCnt(validIdx), normalizedFI(validIdx), 1, color, 'filled');
        s.MarkerFaceAlpha = params.markerAlpha;
        s.MarkerEdgeAlpha = params.markerAlpha;
    end
    
    % Error bar plot
    if size(plotFI, 2) > 0
        validCols = plotFI(2, :) > 0;
        errorbar(exp(plotFI(1, validCols)), plotFI(2, validCols), ...
            plotFI(3, validCols), 'Color', color);
    end
    
    hold off;
    
    % Set properties
    set(gca, 'XScale', 'log');
    ylim(params.yLim);
    xlim([1, params.nRange]);
    set(gca, 'XTick', params.xTicks);
    set(gca, 'XTickLabel', params.xTickLabels);
    xlabel('number of spikes');
    
    if showLabels
        addAxisLabels();
    end
end

function plotAllGroups(plotFI, colors, params)
    % Plot all normalization groups together
    
    hold on;
    for idNormGroup = 1:5
        if size(plotFI{idNormGroup}, 2) > 0
            validCols = plotFI{idNormGroup}(2, :) > 0;
            errorbar(exp(plotFI{idNormGroup}(1, validCols)), ...
                plotFI{idNormGroup}(2, validCols), ...
                plotFI{idNormGroup}(3, validCols), ...
                'Color', colors(idNormGroup, :));
        end
    end
    hold off;
    
    set(gca, 'XScale', 'log');
    ylim(params.yLimAll);
    xlim([1, params.nRange]);
    set(gca, 'XTick', params.xTicks);
    set(gca, 'XTickLabel', params.xTickLabels);
end

function plotAllGroupsOri(plotFI, colors, params)
    % Plot all orientation groups together
    
    hold on;
    for idNormGroup = 1:5
        if size(plotFI{idNormGroup}, 2) > 0
            validCols = plotFI{idNormGroup}(2, :) > 0;
            errorbar(exp(plotFI{idNormGroup}(1, validCols)), ...
                plotFI{idNormGroup}(2, validCols), ...
                plotFI{idNormGroup}(3, validCols), ...
                'Color', colors(idNormGroup, :));
        end
    end
    hold off;
    
    set(gca, 'XScale', 'log');
    ylim(params.yLimAll);
    xlim([1, params.nRange]);
    set(gca, 'XTick', params.xTicks);
    set(gca, 'XTickLabel', params.xTickLabels);
    xlabel('number of spikes');
end

function addAxisLabels()
    % Add y-axis labels
    
    text('Units', 'Normalized', 'Rotation', 90, 'Position', [-0.35, 0.5], ...
        'String', 'Fisher information', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', 'k');
    text('Units', 'Normalized', 'Rotation', 90, 'Position', [-0.25, 0.5], ...
        'String', 'per spike', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', 'k');
end

function addNormIndexLegend(labels, colors)
    % Add legend for normalization index groups
    
    text('Units', 'Normalized', 'Position', [0.8, 0.9], ...
        'String', 'norm. ind.', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', 'k');
    
    for ii = 1:5
        text('Units', 'Normalized', 'Position', [0.8, 0.9 - 0.1*ii], ...
            'String', labels{ii}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, ...
            'Color', colors(ii, :));
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
end
