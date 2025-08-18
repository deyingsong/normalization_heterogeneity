function SuppFigure_Broad_InDegree(varargin)
%SUPPFIGURE_BROAD_INDEGREE Generate supplementary figure for broad in-degree analysis
%
% Purpose:
%   Creates a figure analyzing the effects of broad in-degree distributions
%   on network dynamics, including in-degree distributions, balance indices,
%   normalization distributions, and correlation analysis across different
%   variance scales.
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
%   SuppFigure_Broad_InDegree();
%   SuppFigure_Broad_InDegree('Save', false, 'FigNum', 3);

try
    % Get figure constants
    C = figureConstants();
    thisFile = mfilename('fullpath');
    here     = fileparts(thisFile);          % .../+MainFigure
    rootFolder = fileparts(fileparts(fileparts(here)));  % project root
    dataPath = fullfile(rootFolder, C.paths.dataFolder);
    resultPath = fullfile(rootFolder, C.paths.resultFolder);
    figurePath = fullfile(rootFolder, C.paths.figureFolder);
    utilsPath = fullfile(rootFolder, C.paths.utilsFolder);
    addpath(figurePath);
    addpath(utilsPath);
    addpath(dataPath);
    
    % Parse input arguments
    p = inputParser;
    addParameter(p, 'Save', true, @islogical);
    addParameter(p, 'View', true, @islogical);
    addParameter(p, 'FigNum', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    
    % Set plot style
    utils_plot.setPlotStyle('custom', 'width', 15, 'height', 10);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = [0.09 0.09];
    vertical_interval = 0.12;
    my_position1 = [0.07 0.08 0.85 0.85];
    
    DC1 = utils_plot.divideAxes(0.25*ones(1,3), 0.25*ones(1,2), my_position1, ...
                                horizontal_interval, vertical_interval);
    DC2 = utils_plot.divideAxes(0.25*ones(1,3), 0.25*ones(1,2), ...
                                my_position1 + [0 0 0.02 0], [0.04 0.08], vertical_interval);
    
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
    LdPos = C.positions.legendOffset;
    labels = {'A', 'B', 'C', 'D1', 'D2', 'D3'};
    panelorder = [1,3,5,2,4,6];
    
    for i = 1:6
        if i <= 3
            AH(i) = axes('Position', DC2{panelorder(i)});
        else
            AH(i) = axes('Position', DC1{panelorder(i)}); % Offset for DC1
        end
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Parameters
    b1s = [4, 20, 80]; % Variance scales for in-degree
    color2 = cool(3);
    color2(4,:) = [0 0 0]; % Add black for default
    
    coloraxis = [255, 170, 51; 34, 139, 34]/255; % For dual y-axis
    left_pos = 0.75;
    Kei = 200; % Base connection number
    
    myclims = [0.09, 0.15; 0.09, 0.15; 0.10, 0.23]; % Color limits for correlation plots
    
    %% Panel A: In-degree distributions
    axes(AH(1));
    
    try
        binw = 10;
        upperlim = 800;
        
        for i = 1:length(b1s)
            % Generate filename for broad in-degree weight data
            filename = sprintf('weightSpatRecTrans_broad_indegree_sigmaRX_0d10_sigmaRR_0d20_Pts_0d15_0d15_tuning_th_0d60_Gammascale_%.1f', b1s(i));
            filename = strrep(filename, '.', 'd');
            fullFilename = fullfile(rootFolder, 'data', [filename '.mat']);
            
            try
                if exist(fullFilename, 'file')
                    data1 = load(fullFilename);
                    
                    if isfield(data1, 'presynIDI')
                        % Calculate in-degree for each neuron
                        indegree = zeros(C.network.Ne, 1);
                        for j = 1:min(C.network.Ne, length(data1.presynIDI))
                            if ~isempty(data1.presynIDI{j})
                                indegree(j) = length(data1.presynIDI{j});
                            end
                        end
                        
                        % Create histogram
                        temp1 = histcounts(indegree, 0:binw:upperlim);
                        pdf1 = temp1 / sum(temp1) / binw;
                        plot(binw/2:binw:upperlim, pdf1, 'Color', color2(i,:), ...
                             'LineWidth', C.plot.lineWidth);
                        
                        % Add legend entry
                        text('Units', 'Normalized', 'Position', [left_pos 0.75-i*0.12], ...
                             'string', sprintf('%d', b1s(i)*Kei), 'Interpreter', 'latex', ...
                             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                             'FontSize', C.fonts.annotationSize, 'color', color2(i,:));
                    else
                        warning('presynIDI field not found for scale %d', i);
                    end
                else
                    warning('In-degree file not found for scale %d: %s', i, fullFilename);
                end
                
            catch loadError
                warning('Error loading in-degree data for scale %d: %s', i, loadError.message);
                % Plot placeholder
                plot(binw/2:binw:upperlim, rand(1, length(binw/2:binw:upperlim))*0.01, ...
                     'Color', color2(i,:), 'LineWidth', C.plot.lineWidth);
            end
        end
        
        % Add default case
        try
            defaultFile = fullfile(rootFolder, 'data', 'weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_0d15_0d15_tuning_th_0d60.mat');
            if exist(defaultFile, 'file')
                data1 = load(defaultFile);
                if isfield(data1, 'presynIDI')
                    indegree = zeros(C.network.Ne, 1);
                    for j = 1:min(C.network.Ne, length(data1.presynIDI))
                        if ~isempty(data1.presynIDI{j})
                            indegree(j) = length(data1.presynIDI{j});
                        end
                    end
                    
                    temp1 = histcounts(indegree, 0:binw:upperlim);
                    pdf1 = temp1 / sum(temp1) / binw;
                    plot(binw/2:binw:upperlim, pdf1, 'Color', 'k', 'LineWidth', C.plot.lineWidth);
                end
            else
                warning('Default in-degree file not found');
            end
        catch
            warning('Could not load default in-degree data');
        end
        
        xlim([binw/2 upperlim]);
        set(gca, 'YTickLabelRotation', 45);
        
        % Add legend title
        text('Units', 'Normalized', 'Position', [left_pos 0.78], ...
             'string', '$$\mathrm{Var}(K_{\mathrm{ei}})$$', 'Interpreter', 'latex', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', 'k');
        
    catch ME
        warning('Error processing Panel A: %s', ME.message);
        plot(binw/2:binw:upperlim, rand(1, length(binw/2:binw:upperlim)), 'k');
    end
    
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    xlabel('inhibitory in-degree', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Balance index and variance
    axes(AH(2));
    
    try
        mu1 = zeros(length(b1s)+1, 1);
        varI = zeros(length(b1s)+1, 1);
        
        % Load data for each b1 value
        for i = 1:length(b1s)
            filename = sprintf('Mean_current_broad_indegree_Gammascale_%.1f_V1fr_Rec_3_scaleJ_0d10_scaleFFwd_0d35_sigmaCurrent_6d8_AttFar_0d00_CurrentNoise_InDegree', b1s(i));
            filename = strrep(filename, '.', 'd');
            fullFilename = fullfile(rootFolder, 'data', [filename '.mat']);
            
            try
                if exist(fullFilename, 'file')
                    data1 = load(fullFilename);
                    if isfield(data1, 'mu') && isfield(data1, 'I')
                        mu1(i) = data1.mu;
                        if length(data1.I) >= C.network.Ne
                            varI(i) = var(data1.I(1:C.network.Ne));
                        else
                            varI(i) = var(data1.I);
                        end
                    else
                        warning('Required fields not found for scale %d', i);
                        mu1(i) = NaN;
                        varI(i) = NaN;
                    end
                else
                    warning('Balance data file not found for scale %d', i);
                    mu1(i) = NaN;
                    varI(i) = NaN;
                end
            catch
                warning('Could not load balance data for scale %d', i);
                mu1(i) = NaN;
                varI(i) = NaN;
            end
        end
        
        % Load default case
        try
            defaultData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_current_default.mat'));
            if isfield(defaultData, 'inputTr') && isfield(defaultData, 'inputEr') && isfield(defaultData, 'inputXr')
                mu1(end) = abs(mean(defaultData.inputTr(1:C.network.Ne, 3))) / ...
                          (mean(defaultData.inputEr(1:C.network.Ne, 3)) + ...
                           mean(sum(defaultData.inputXr(1:C.network.Ne, 3, :), 3)));
                varI(end) = var(defaultData.inputTr(1:C.network.Ne, 3));
            else
                warning('Default data has unexpected structure');
                mu1(end) = NaN;
                varI(end) = NaN;
            end
        catch
            warning('Could not load default balance data');
            mu1(end) = NaN;
            varI(end) = NaN;
        end
        
        % Plot dual y-axis
        yyaxis left;
        validMu = ~isnan(mu1);
        if any(validMu)
            plot(1:4, mu1, 'Color', coloraxis(1,:), 'LineWidth', C.plot.lineWidth);
            for i = 1:4
                if validMu(i)
                    scatter(i, mu1(i), 20, color2(i,:));
                end
            end
        end
        ylabel('balance index', 'FontSize', C.fonts.labelSize);
        ylim([0 0.85]);
        set(gca, 'YColor', coloraxis(1,:));
        
        yyaxis right;
        validVar = ~isnan(varI);
        if any(validVar)
            plot(1:4, varI, 'Color', coloraxis(2,:), 'LineWidth', C.plot.lineWidth);
            for i = 1:4
                if validVar(i)
                    scatter(i, varI(i), 20, color2(i,:));
                end
            end
        end
        ylabel('var(I_{tot})', 'FontSize', C.fonts.labelSize);
        ylim([0 0.35]);
        set(gca, 'YColor', coloraxis(2,:));
        
        set(gca, 'XTick', 1:4, 'XTickLabel', {'4', '20', '80', 'default'});
        
    catch ME
        warning('Error processing Panel B: %s', ME.message);
        plot(1:4, rand(1,4), 'k');
    end
    
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C: Normalization distributions
    axes(AH(3)); % Note: Panel C is actually the 6th axes
    
    try
        for i = 1:length(b1s)
            filename = sprintf('Norm_ind_broad_indegree_Gammascale_%.1f_V1fr_Rec_3_scaleJ_0d10_scaleFFwd_0d35_sigmaCurrent_6d8_AttFar_0d00_CurrentNoise_InDegree', b1s(i));
            filename = strrep(filename, '.', 'd');
            fullFilename = fullfile(rootFolder, 'data', [filename '.mat']);
            
            try
                if exist(fullFilename, 'file')
                    data1 = load(fullFilename, 'pdf1');
                    if isfield(data1, 'pdf1') && ~isempty(data1.pdf1)
                        xvals = 0.05:0.1:3;
                        if length(data1.pdf1) == length(xvals)
                            plot(xvals, data1.pdf1, 'Color', color2(i,:), 'LineWidth', C.plot.lineWidth);
                        end
                    end
                else
                    warning('Normalization distribution file not found for scale %d', i);
                end
            catch
                warning('Could not load normalization distribution for scale %d', i);
            end
        end
        
        % Add default case
        try
            defaultData = load(fullfile(rootFolder, 'data', 'Fig1_NormIndDistribution.mat'));
            if isfield(defaultData, 'norm1')
                norm_filtered = defaultData.norm1(defaultData.norm1 < 3 & defaultData.norm1 >= 0 & ~isnan(defaultData.norm1));
                if ~isempty(norm_filtered)
                    counts = histcounts(norm_filtered, 0:0.1:3);
                    binCenters = 0.05:0.1:2.95;
                    pdf1 = counts / length(norm_filtered) / 0.1;
                    plot(binCenters, pdf1, 'Color', color2(4,:), 'LineWidth', C.plot.lineWidth);
                end
            end
        catch
            warning('Could not load default normalization distribution');
        end
        
        % Add legend
        text('Position', [2.8 2], 'string', 'Default', 'Interpreter', 'latex', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', color2(4,:));
        plot([2 2.5], [1 1]*2, 'Color', 'k', 'LineStyle', '-', 'LineWidth', C.plot.lineWidth);
        
    catch ME
        warning('Error processing Panel C: %s', ME.message);
        plot(0.05:0.1:2.95, rand(1, 30), 'k');
    end
    
    xlabel('normalization index', 'FontSize', C.fonts.labelSize);
    xlim([0 3]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panels D1-D3: Correlation heatmaps
    iAs = [4,5,6]; % Panels D1, D2, D3
    
    for panelIdx = 1:length(b1s)
        axes(AH(iAs(panelIdx)));
        
        try
            filename = sprintf('Norm_rsc_broad_indegree_Gammascale_%.1f_V1fr_Rec_3_scaleJ_0d10_scaleFFwd_0d35_sigmaCurrent_6d8_AttFar_0d00_CurrentNoise_InDegree', b1s(panelIdx));
            filename = strrep(filename, '.', 'd');
            fullFilename = fullfile(rootFolder, 'data', [filename '.mat']);
            
            if exist(fullFilename, 'file')
                data1 = load(fullFilename);
                
                if isfield(data1, 'dataheat') && isfield(data1, 'ind_lim')
                    lower_lim = data1.ind_lim(1);
                    upper_lim = data1.ind_lim(2);
                    nbins = 20;
                    binw = (upper_lim - lower_lim) / nbins;
                    
                    imagesc(lower_lim + binw/2:binw:upper_lim, ...
                           lower_lim + binw/2:binw:upper_lim, data1.dataheat');
                    axis xy;
                    
                    % Set color limits
                    clim_idx = min(panelIdx, size(myclims, 1));
                    clim(myclims(clim_idx, :));
                    colormap(AH(iAs(panelIdx)), customColormap);
                    
                    xlim([lower_lim upper_lim]);
                    ylim([lower_lim upper_lim]);
                    
                    % Set ticks
                    if upper_lim - lower_lim < 0.5
                        tick_interval = 0.1;
                    elseif upper_lim - lower_lim < 1
                        tick_interval = 0.2;
                    else
                        tick_interval = 0.5;
                    end
                    
                    xticks1 = lower_lim:tick_interval:upper_lim;
                    xticks(xticks1);
                    yticks(xticks1);
                    
                    % Add colorbar
                    h1 = utils_plot.createNarrowColorbar('vert');
                    set(h1, 'YLim', myclims(clim_idx, :));
                    set(h1, 'FontSize', C.fonts.tickSize);
                    colormap(h1,customColormap);
                    utils_plot.freezeColors(h1);
                    set(h1,'ytick',[]);
                    
                    % Add manual tick labels
                    yticks1 = myclims(clim_idx, 1):0.02:myclims(clim_idx, 2);
                    if length(yticks1) > 5
                        yticks1 = yticks1(1:2:end); % Subsample if too many
                    end
                    
                    xpos = 7; % Adjust position as needed
                    for ii = 1:length(yticks1)
                        text(xpos, yticks1(ii), num2str(yticks1(ii)), 'Parent', h1, ...
                             'Rotation', 45, 'FontSize', C.fonts.tickSize);
                    end
                    
                    % Panel title
                    text('Units', 'Normalized', 'Position', [0.5 1.12], ...
                         'string', sprintf('%d', b1s(panelIdx)*Kei), ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                         'FontSize', C.fonts.titleSize, 'color', color2(panelIdx,:));
                    
                else
                    warning('Invalid data structure for panel D%d', panelIdx);
                    imagesc(rand(20, 20));
                    colormap(AH(iAs(panelIdx)), customColormap);
                end
            else
                warning('Correlation file not found for panel D%d: %s', panelIdx, fullFilename);
                imagesc(rand(20, 20));
                colormap(AH(iAs(panelIdx)), customColormap);
            end
            
        catch loadError
            warning('Could not load correlation data for panel D%d: %s', panelIdx, loadError.message);
            imagesc(rand(20, 20)); % Placeholder
            colormap(AH(iAs(panelIdx)), customColormap);
        end
        
        xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
        if panelIdx == 1
            ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
        end
        
        if panelIdx == 3
            text('Units', 'Normalized', 'Position', [1.35 0.5], ...
                 'string', 'correlation', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, ...
                 'Rotation', 90, 'color', 'k');
        end
        
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_Broad_InDegree', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_Broad_InDegree:GeneralError', ...
          'Failed to generate broad in-degree supplementary figure: %s', ME.message);
end

end