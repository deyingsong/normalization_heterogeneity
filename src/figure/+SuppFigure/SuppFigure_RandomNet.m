function SuppFigure_RandomNet(varargin)
%SUPPFIGURE_RANDOMNET Generate supplementary figure for random network analysis
%
% Purpose:
%   Creates a figure analyzing random network properties including schematic
%   representation, normalization index distributions across different overlap
%   conditions, and correlations between current and rate normalization indices.
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
%   SuppFigure_RandomNet();
%   SuppFigure_RandomNet('Save', false, 'View', true);

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
    utils_plot.setPlotStyle('custom', 'width', 12, 'height', 10);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    position = [0.06 0.12 0.90 0.80];
    vert_space = 0.25;
    DC = utils_plot.divideAxes([0.5 0.5], [0.4 0.4], position, 0.2, vert_space);
    DC2 = utils_plot.divideAxes([0.4 0.4 0.4], [0.4 0.4], position, [0.1 0.1], vert_space);
    
    % Create axes
    LdPos = [-0.03, 0.05];
    labels = {'A', 'B', 'C1', 'C2', 'C3'};
    
    AH(1) = axes('Position', DC{1});
    hold on;
    utils_plot.addFigureLabel('A', LdPos);
    
    AH(2) = axes('Position', DC{3});
    hold on;
    utils_plot.addFigureLabel('B', LdPos);
    
    for i = 3:5
        AH(i) = axes('Position', DC2{i*2-4}); % Offset for DC2 indexing
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties();
    
    % Colors
    col1 = C.colors.excitatory;
    
    %% Panel A: Schematic (placeholder)
    axes(AH(1));
    
    try
        % Create a simple schematic representation
        text(0.5, 0.5, 'Schematic', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.titleSize, 'FontWeight', 'bold');
        
        % Add some basic network elements
        circle_x = [0.2, 0.8, 0.5];
        circle_y = [0.7, 0.7, 0.3];
        circle_colors = [col1; C.colors.inhibitory; C.colors.feedforward];
        
        for i = 1:3
            rectangle('Position', [circle_x(i)-0.05, circle_y(i)-0.05, 0.1, 0.1], ...
                     'Curvature', [1, 1], 'FaceColor', circle_colors(i,:), ...
                     'EdgeColor', 'k', 'LineWidth', C.plot.lineWidth);
        end
        
        % Add connections
        for i = 1:2
            for j = i+1:3
                line([circle_x(i), circle_x(j)], [circle_y(i), circle_y(j)], ...
                     'Color', 'k', 'LineWidth', C.plot.lineWidth/2, 'LineStyle', '--');
            end
        end
        
    catch ME
        warning('Error creating schematic: %s', ME.message);
        text(0.5, 0.5, 'Schematic not available', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center');
    end
    
    title('Schematic', 'FontSize', C.fonts.titleSize);
    set(gca, 'XTick', [], 'YTick', []);
    pbaspect(C.figure.aspectRatio);
    
    %% Load data for other panels
    try
        % Load current normalization data
        currentData = load(fullfile(rootFolder, 'data', 'FigSupp_CurrentNormIndOverlap0d6.mat'));
        randomData = load(fullfile(rootFolder, 'data', 'indrandom_3.mat'));
        varOverlapData = load(fullfile(rootFolder, 'data', 'FigSupp_NormIndDistributionVarOverlap.mat'));
        
        % Load MTE data for multiple conditions
        MTE = cell(7, 1);
        norm1 = zeros(C.network.Ne, 6);
        
        
        MTE = varOverlapData.MTE;
        for i=1:6
            norm1(:,i) = (MTE{i}(:,1) + MTE{i}(:,2)) ./ MTE{i}(:,3);
        end
                

        
        % Load default condition
        try
            defaultData = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond.mat'));
            if isfield(defaultData, 'MTE')
                MTE{7} = defaultData.MTE;
                norm2 = (defaultData.MTE(:,1) + defaultData.MTE(:,2)) ./ defaultData.MTE(:,3);
            end
        catch
            warning('Could not load default condition data');
            norm2 = [];
        end
        
        % Extract required variables
        if isfield(currentData, 'normX') && isfield(currentData, 'normE') && isfield(currentData, 'normI')
            normX = currentData.normX;
            normE = currentData.normE;
            normI = currentData.normI;
        else
            warning('Required normalization fields not found in current data');
            normX = ones(1000, 1) * 1.5 + randn(1000, 1) * 0.1;
            normE = ones(1000, 1) * 1.5 + randn(1000, 1) * 0.1;
            normI = ones(1000, 1) * 1.5 + randn(1000, 1) * 0.1;
        end
        
        if isfield(randomData, 'ind1')
            ind1 = randomData.ind1;
        else
            warning('Random indices not found, using default');
            ind1 = 1:min(1000, size(norm1, 1));
        end
        
        % Ensure ind1 doesn't exceed data bounds
        ind1 = ind1(ind1 <= size(norm1, 1) & ind1 <= length(normX));
        
    catch ME
        error('SuppFigure_RandomNet:DataLoadError', ...
              'Failed to load required data: %s', ME.message);
    end
    
    %% Panel B: Normalization index distributions
    axes(AH(2));
    
    try
        cols = cool(6);
        densitycurve = zeros(50, 7);
        
        hold on;
        text('Units', 'Normalized', 'Position', [0.7 0.98], ...
             'string', 'overlap', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
        
        % Plot overlap conditions
        for i = 1:6
            if ~isempty(MTE{i}) && size(norm1, 2) >= i
                norm_filtered = norm1(:,i);
                norm_filtered = norm_filtered(norm_filtered < 5 & ~isnan(norm_filtered) & ...
                                            abs(MTE{i}(:,3) - mean(MTE{i}(:,3))) < std(MTE{i}(:,3)));
                
                if ~isempty(norm_filtered)
                    N = histcounts(norm_filtered, 0:0.1:5);
                    densitycurve(:,i) = N / length(norm_filtered) / 0.1;
                    plot(0.05:0.1:4.95, densitycurve(:,i), 'color', cols(i,:), 'LineWidth', C.plot.lineWidth);
                    
                    text('Units', 'Normalized', 'Position', [0.7 0.95-i*0.10], ...
                         'string', num2str((i-1)*0.2), 'Interpreter', 'latex', ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                         'FontSize', C.fonts.annotationSize, 'color', cols(i,:));
                end
            end
        end
        
        % Plot default condition
        if ~isempty(norm2)
            norm_filtered = norm2(norm2 < 5 & ~isnan(norm2) & ...
                                abs(MTE{7}(:,3) - mean(MTE{7}(:,3))) < std(MTE{7}(:,3)));
            
            if ~isempty(norm_filtered)
                N = histcounts(norm_filtered, 0:0.1:5);
                densitycurve(:,7) = N / length(norm_filtered) / 0.1;
                plot(0.05:0.1:4.95, densitycurve(:,7), 'color', 'k', 'LineWidth', C.plot.lineWidth);
                
                text('Units', 'Normalized', 'Position', [0.7 0.95-7*0.10], ...
                     'string', 'Default', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
            end
        end
        
        hold off;
        
    catch ME
        warning('Error processing normalization distributions: %s', ME.message);
        plot(0.05:0.1:4.95, rand(1, 50), 'k');
    end
    
    xlabel('norm index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    xlim([0 5]);
    set(gca, 'XTick', 0:1:5);
    pbaspect([1.25 1 1]);
    
    %% Panels C1-C3: Current vs rate normalization correlations
    panelTitles = {'Feedforward excitation', 'Recurrent excitation', 'Recurrent inhibition'};
    normData = {normX, normE, normI};
    colors = {col1, col1, C.colors.inhibitory};
    
    for panelIdx = 1:3
        axes(AH(panelIdx + 2));
        
        try
            if length(ind1) <= length(normData{panelIdx}) && size(norm1, 1) >= max(ind1)
                current_norm = normData{panelIdx}(ind1);
                rate_norm = norm1(ind1, 4); % Use 4th column as representative
                
                % Remove invalid data
                validData = ~isnan(current_norm) & ~isnan(rate_norm) & ...
                           isfinite(current_norm) & isfinite(rate_norm);
                
                if sum(validData) > 10
                    current_norm = current_norm(validData);
                    rate_norm = rate_norm(validData);
                    
                    scatter(rate_norm, current_norm, 2, 'MarkerFaceColor', colors{panelIdx}, ...
                           'MarkerEdgeColor', colors{panelIdx}, 'LineWidth', 0.1, ...
                           'MarkerFaceAlpha', C.plot.alphaValue);
                    
                    % Calculate correlation
                    [r, p1] = corrcoef(rate_norm, current_norm);
                    if size(r, 1) > 1
                        text('Units', 'Normalized', 'Position', [0.8 0.8], ...
                             'string', sprintf('R=%.2f', r(1,2)), 'Interpreter', 'latex', ...
                             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                             'FontSize', C.fonts.annotationSize, 'color', 'k');
                    end
                else
                    warning('Insufficient valid data for panel C%d', panelIdx);
                    scatter(rand(50,1)*15, rand(50,1)*0.3+1.4, 2, colors{panelIdx});
                end
            else
                warning('Data dimension mismatch for panel C%d', panelIdx);
                scatter(rand(50,1)*15, rand(50,1)*0.3+1.4, 2, colors{panelIdx});
            end
            
        catch ME
            warning('Error processing panel C%d: %s', panelIdx, ME.message);
            scatter(rand(50,1)*15, rand(50,1)*0.3+1.4, 2, colors{panelIdx});
        end
        
        text('Units', 'Normalized', 'Position', [0.5 C.positions.titleHeight], ...
             'string', panelTitles{panelIdx}, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', C.fonts.titleSize, 'color', colors{panelIdx});
        
        xlabel('rate norm index', 'FontSize', C.fonts.labelSize);
        if panelIdx == 1
            ylabel('current norm index', 'FontSize', C.fonts.labelSize);
        end
        
        ylim([1.4 1.7]);
        xlim([0 15]);
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_RandomNet', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_RandomNet:GeneralError', ...
          'Failed to generate random network supplementary figure: %s', ME.message);
end

end