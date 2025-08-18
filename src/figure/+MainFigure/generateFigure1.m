function generateFigure1()
% GENERATEFIGURE1_REVISED Generate Figure 1 for the normalization paper
%
% Purpose:
%   Creates a multi-panel figure showing:
%   - Orientation preference map with connectivity
%   - Gabor image demonstrations
%   - Firing rate patterns
%   - Normalization index distributions
%
% Inputs:
%   None (loads data from files)
%
% Outputs:
%   Saves figure as PDF in results folder
%
% Usage Example:
%   generateFigure1_revised();

    % Initialize
    try
        % Load constants
        C = figureConstants();
        
        % Setup paths
        dataPath = fullfile(pwd, C.paths.dataFolder);
        resultPath = fullfile(pwd, C.paths.resultFolder);
        
        % Create result directory if it doesn't exist
        if ~exist(resultPath, 'dir')
            mkdir(resultPath);
        end
        
        % Setup plot options
        utils_plot.setPlotStyle('custom', 'path', resultPath, ...
            'width', 16.4, 'height', 8);
        
        % Prepare figure
        hFig = figure(1);
        clf;
        utils_plot.matchFigureAspectRatio();
        
        % Create panel layout
        DC = utils_plot.divideAxes([0.45, 0.18, 0.45, 0.4], ...
            [0.3, 0.3, 0.3], [0.02, 0.1, 0.82, 0.82], ...
            [0.14, 0.01, 0.03], [0.05, 0.05])';
        DC1 = utils_plot.divideAxes([0.40, 0.75, 0.4], ...
            [0.45, 0.45], [0.12, 0.13, 0.82, 0.75], ...
            [0.05, 0.05], 0.25)';
        
        % Rearrange panels
        DCf = arrangePanels(DC, DC1);
        
        % Create axes
        AH = createAxes(DCf);
        utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
        % Plot each panel
        plotPanelA(AH(1), AH(2), dataPath, C, DCf);  % Orientation map
        plotPanelsB(AH(3:5), dataPath, C);       % Gabor images
        plotFiringRates(AH(6:8), dataPath, C, DCf);   % Firing rate patterns
        plotPanelC(AH(9), dataPath, C);          % Model norm index
        plotPanelD(AH(10), dataPath, C);         % Experimental data
        
        % Save figure
        saveFigure(hFig, resultPath, 'Fig1');
        
    catch ME
        fprintf('Error in generateFigure1_revised: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s, Line: %d\n', ME.stack(i).file, ME.stack(i).line);
        end
        rethrow(ME);
    end
end

function DCf = arrangePanels(DC, DC1)
    % Arrange panel positions
    DCf = cell(10, 1);
    DCf{1} = DC{1};
    DCf{2} = DC{9};
    DCf{3} = DC{2};
    DCf{4} = DC{6};
    DCf{5} = DC{10};
    DCf{6} = DC{3};
    DCf{7} = DC{7};
    DCf{8} = DC{11};
    DCf{9} = DC1{3};
    DCf{10} = DC1{6};
end

function AH = createAxes(DCf)
    % Create axes for all panels
    nPanels = length(DCf);
    AH = gobjects(nPanels, 1);
    
    labels = {'A', '', 'B1', 'B2', 'B3', '', '', '', 'C', 'D'};
    labelPos = [-0.03, 0.04];
    labelPos1 = [0, 0.20];
    
    for i = 1:nPanels
        AH(i) = axes('Position', DCf{i});
        hold(AH(i), 'on');
        
        % Add labels where appropriate
        if ~isempty(labels{i})
            if i == 1
                utils_plot.addFigureLabel(labels{i}, [0.01, 0]);
            elseif i >= 3 && i <= 5
                utils_plot.addFigureLabel(labels{i}, labelPos1);
            elseif i >= 9
                utils_plot.addFigureLabel(labels{i}, labelPos);
            end
        end
        set(0, 'DefaultAxesTitleFontWeight', 'normal');
    end
end

function plotPanelA(ax1, ax2, dataPath, C, DCf)
    % Plot orientation preference map with connectivity
    
    axes(ax1);
    axis off;
    
    axes(ax2);
    
    % Load data with error handling
    try
        data1 = load(fullfile(dataPath, 'theta_map_MTE_two_rec.mat'));
        data2 = load(fullfile(dataPath, ...
            'weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_0d15_0d15_tuning_th_0d60.mat'));
    catch ME
        error('Failed to load orientation map data: %s', ME.message);
    end
    
    % Validate data
    if ~isfield(data1, 'theta_mapE')
        error('theta_mapE not found in data file');
    end
    if ~isfield(data2, 'Wrr')
        error('Wrr not found in weight data file');
    end
    
    % Plot orientation map
    imagesc(data1.theta_mapE);
    colormap(ax2, 'hsv');
    
    % Set dimensions
    w_h = 0.115 * [2, 16.4/8];
    set(gca, 'Units', 'Normalized', 'Position', [0.015, 0.1, w_h]);
    set(gca, 'DataAspectRatio', C.figure.aspectRatio);
    
    % Add connectivity overlay
    pre_ID = 10110;
    Ney = 200;
    Kee_eets = 200;
    Ke = 350;
    
    if pre_ID > 0 && pre_ID * Ke <= length(data2.Wrr)
        post_ID = double(data2.Wrr((pre_ID-1)*Ke+1:(pre_ID-1)*Ke+Kee_eets));
        scatter(mod(post_ID-1, Ney)+1, ceil(post_ID/Ney), 3, 'k', 'filled');
        scatter(mod(pre_ID-1, Ney)+1, ceil(pre_ID/Ney), 9, 'w', 'filled');
    else
        warning('Invalid pre_ID or weight matrix size');
    end
    
    axis off;
    xlim([0, 200]);
    ylim([0, 100]);
    
    % Add colorbar
    h1 = colorbar;
    set(h1, 'Position', [0.25, 0.1, 0.005, w_h(2)]);
    caxis([0, 1]);
    
    % Add colorbar label
    text('Units', 'Normalized', 'Position', [1.2, 0.5], ...
        'String', 'Preferred ori. (deg.)', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', C.fonts.annotationSize, ...
        'Rotation', 90);
    
    % Add colorbar ticks
    set(h1, 'YTick', [0, 0.5, 1], 'YTickLabel', {}, 'FontSize', C.fonts.tickSize);
    
    % Add tick labels manually
    yticks1 = [0, 0.5, 1];
    yticklabels1 = {'0', '90', '180'};
    for i = 1:length(yticks1)
        text('Units', 'Normalized', ...
            'Position', [1.10, 0.03+(i-1)*w_h(2)/2/DCf{2}(4)], ...
            'String', yticklabels1{i}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', C.fonts.tickSize);
    end
end

function plotPanelsB(axes_handles, dataPath, C)
    % Plot Gabor image demonstrations
    
    try
        data = load(fullfile(dataPath, 'Fig1_GaborImageDemonstration.mat'));
    catch ME
        error('Failed to load Gabor image data: %s', ME.message);
    end
    
    if ~isfield(data, 'image1') || length(data.image1) < 3
        error('Invalid Gabor image data structure');
    end
    
    legendPos = [0.8, 1.6];
    labels = {'$$r_1$$', '$$r_2$$', '$$r_{1+2}$$'};
    
    for i = 1:3
        axes(axes_handles(i));
        
        % Process image - shift quadrants for display
        if size(data.image1{i}, 1) >= 58 && size(data.image1{i}, 2) >= 29
            temp = zeros(29, 58);
            temp(1:29, 1:29) = data.image1{i}(30:58, 1:29);
            temp(1:29, 30:58) = data.image1{i}(1:29, 1:29);
            
            imagesc(temp);
            colormap(axes_handles(i), 'gray');
            axis off;
            box off;
            set(gca, 'DataAspectRatio', C.figure.aspectRatio);
            set(gca, 'YDir', 'reverse');
            
            % Add label
            text('Units', 'Normalized', 'Position', legendPos, ...
                'String', labels{i}, ...
                'Interpreter', 'latex', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', C.fonts.annotationSize);
        else
            warning('Image %d has unexpected dimensions', i);
        end
    end
end

function plotFiringRates(axes_handles, dataPath, C, DCf)
    % Plot firing rate patterns
    
    try
        data = load(fullfile(dataPath, 'Fig1_Pattern3cond.mat'));
    catch ME
        error('Failed to load firing rate pattern data: %s', ME.message);
    end
    
    if ~isfield(data, 'MTE') || size(data.MTE, 1) ~= 20000
        error('Invalid MTE data structure');
    end
    
    w_h = 0.115 * [2, 16.4/8];
    left_pos = 0.4053;
    height_pos = [0.674, 0.387, 0.1];
    
    for i = 1:3
        axes(axes_handles(i));
        
        % Reshape and plot firing rate data
        frData = reshape(data.MTE(:, i), 200, 100)';
        imagesc(frData);
        
        % Set colormap
        oldcmap = colormap(axes_handles(i), 'gray');
        colormap(axes_handles(i), flipud(oldcmap));
        
        % Set position and properties
        set(gca, 'Units', 'Normalized', ...
            'Position', [left_pos, height_pos(i), w_h]);
        set(gca, 'DataAspectRatio', C.figure.aspectRatio);
        caxis(C.analysis.firingRateRange);
        
        axis off;
        xlim([0, 200]);
        ylim([0, 100]);
    end
    
    % Add shared colorbar
    h2 = colorbar;
    h2_height = 0.81;
    set(h2, 'Position', [0.643, 0.1, 0.008, h2_height]);
    caxis(C.analysis.firingRateRange);
    
    % Add colorbar label
    text('Units', 'Normalized', 'Position', [1.21, 1.5], ...
        'String', 'Firing rate (Hz)', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', C.fonts.annotationSize, ...
        'Rotation', 90);
    
    % Add colorbar ticks
    set(h2, 'YTick', 0:10:90, 'YTickLabel', {}, ...
        'FontSize', C.fonts.tickSize);
    
    % Add tick labels manually
    yticks1 = 0:10:90;
    for i = 1:length(yticks1)
        text('Units', 'Normalized', ...
            'Position', [1.12, (i-1)*(h2_height+0.02)/9/DCf{7}(4)], ...
            'String', num2str(yticks1(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', C.fonts.tickSize);
    end
end

function plotPanelC(ax, dataPath, C)
    % Plot model normalization index distribution
    
    axes(ax);
    
    try
        data = load(fullfile(dataPath, 'Fig1_NormIndDistribution.mat'));
    catch ME
        error('Failed to load normalization index data: %s', ME.message);
    end
    
    if ~isfield(data, 'norm1')
        error('norm1 field not found in data');
    end
    
    % Process data
    temp = data.norm1(data.norm1 <= 3);
    if isempty(temp)
        warning('No valid normalization indices found');
        return;
    end
    
    % Calculate histogram
    edges = C.analysis.normIndexBins;
    counts = histcounts(temp, edges);
    binWidth = edges(2) - edges(1);
    density = counts / length(temp) / binWidth;
    binCenters = edges(1:end-1) + binWidth/2;
    
    % Plot
    plot(binCenters, density, 'k', 'LineWidth', C.plot.lineWidth);
    xlim(C.analysis.normIndexRange);
    xlabel('norm index');
    ylabel('probability density');
    title('Model');
    
    % Add equation
    text('Units', 'Normalized', 'Position', [0.85, 0.7], ...
        'String', ['$$\mathrm{norm \ index}$$' newline newline ...
                   '$$=\frac{r_1+r_2}{r_{1+2}}$$'], ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', C.fonts.annotationSize);
end

function plotPanelD(ax, dataPath, C)
    % Plot experimental normalization index distributions
    
    axes(ax);
    hold on;
    
    % Load V4 and MT data
    try
        V4data = load(fullfile(dataPath, 'Fig_ExpData_V4_NormIndDistribution.mat'));
        MTdata = load(fullfile(dataPath, 'Fig_ExpData_MT_NormIndDistribution.mat'));
    catch ME
        error('Failed to load experimental data: %s', ME.message);
    end
    
    % Validate data
    if ~isfield(MTdata, 'norminds') || ~isfield(V4data, 'norminds')
        error('norminds field missing in experimental data');
    end
    
    legendPosY = [1.7, 1.3];
    edges = C.analysis.normIndexBins;
    binWidth = edges(2) - edges(1);
    binCenters = edges(1:end-1) + binWidth/2;
    
    % Plot MT data
    temp = MTdata.norminds(MTdata.norminds <= 3);
    if ~isempty(temp)
        counts = histcounts(temp, edges);
        density = counts / length(temp) / binWidth;
        plot(binCenters, density, 'k:', 'LineWidth', C.plot.lineWidth);
        plot([1.7, 2.4], legendPosY(1) * [1, 1], 'k:', ...
            'LineWidth', C.plot.lineWidth);
        text('Position', [3, legendPosY(1)], 'String', 'MT data', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', C.fonts.annotationSize);
    end
    
    % Plot V4 data
    temp = V4data.norminds(V4data.norminds <= 3);
    if ~isempty(temp)
        counts = histcounts(temp, edges);
        density = counts / length(temp) / binWidth;
        plot(binCenters, density, 'k-.', 'LineWidth', C.plot.lineWidth);
        plot([1.7, 2.4], legendPosY(2) * [1, 1], 'k-.', ...
            'LineWidth', C.plot.lineWidth);
        text('Position', [3, legendPosY(2)], 'String', 'V4 data', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', C.fonts.annotationSize);
    end
    
    xlim(C.analysis.normIndexRange);
    ylim([0, 2]);
    xlabel('norm index');
    ylabel('probability density');
    title('Experimental data');
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
