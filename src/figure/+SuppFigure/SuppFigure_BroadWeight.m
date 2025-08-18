function SuppFigure_BroadWeight(varargin)
%SUPPFIGURE_BROADWEIGHT Generate supplementary figure for broad weight distribution analysis
%
% Purpose:
%   Creates a comprehensive figure showing the effects of broad weight distributions
%   on network dynamics, including weight distributions, balance indices, normalization
%   index distributions, and correlation analysis across different variance scales.
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
%   SuppFigure_BroadWeight();
%   SuppFigure_BroadWeight('Save', false, 'FigNum', 3);

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
    utils_plot.setPlotStyle('custom', 'width', 18, 'height', 10);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = [0.09 0.09 0.09];
    vertical_interval = 0.12;
    my_position1 = [0.06 0.08 0.85 0.85];
    
    DC1 = utils_plot.divideAxes(0.25*ones(1,4), 0.25*ones(1,2), my_position1, ...
                                horizontal_interval, vertical_interval);
    DC2 = utils_plot.divideAxes(0.25*ones(1,4), 0.25*ones(1,2), ...
                                my_position1 + [0 0 0.05 0], [0.1 0.2 0.1], vertical_interval);
    
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
    LdPos = C.positions.legendOffset+[0 0.03];
    labels = {'A', 'B', 'C', 'D1', 'D2', 'D3','D4'};
    
    for i = 1:7
        if i <= 3
            AH(i) = axes('Position', DC2{2*i-1});
        else
            AH(i) = axes('Position', DC1{i*2-6}); % Offset for DC2
        end
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    
    % Parameters
    b1s = [0.2, 0.5, 5.0, 8.3]; % Variance scales
    scaleJ = 1/0.1;
    scaledown_ffwd = 1/0.34;
    sigma_current = 6.8;
    att_cond = 'AttFar';
    sim_flag = 'CurrentNoise_BroadWeight';
    
    color2 = cool(4);
    color2(5,:) = [0 0 0]; % Add black for default
    
    coloraxis = [255, 170, 51; 34, 139, 34]/255; % For dual y-axis
    legend_posX = 0.65;
    Jei = 240/10/sqrt(C.network.N);
    
    myclims = [0.02, 0.08; 0.02, 0.08; 0.00, 0.04; 0, 0.02];
    
    %% Panel A: Weight distribution plots
    axes(AH(1));
    
    try
        upperlim = 400;
        binw = 0.1;
        
        for i = 1:length(b1s)
            b1 = b1s(i);
            
            % Generate filename for weight data
            filename = sprintf('RF2D2layer_wide_V1fr_Rec_3_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%s_bgamma_ex_%.1f_ix_%.1f_ee_%.1f_ie_%.1f_ei_%.1f_ii_%.1f', ...
                              1/scaleJ, 1/scaledown_ffwd, sigma_current, att_cond, sim_flag, ...
                              b1, b1, b1, b1, b1, b1);
            filename = strrep(filename, '.', 'd');
            fullFilename = fullfile(rootFolder, 'data', [filename '.mat']);
            
            try
                if exist(fullFilename, 'file')
                    data1 = load(fullFilename, 'Jrr', 'param');
                    
                    if isfield(data1, 'Jrr') && isfield(data1, 'param')
                        % Extract connection weights
                        m = data1.param.Kr(1,2);
                        n = data1.param.Kr(2,2);
                        K = data1.param.Ni;
                        
                        Xm_idx = reshape(bsxfun(@plus, (0:K-1)' * (m+n), 1:m), [], 1) + ...
                                data1.param.Ne * sum(data1.param.Kr(:,1));
                        
                        weights = data1.Jrr(Xm_idx) / data1.param.Jr(1,2);
                        weights = weights(weights > 0 & weights < upperlim); % Filter valid weights
                        
                        if ~isempty(weights)
                            counts = histcounts(weights, 0:binw:upperlim);
                            binCenters = binw/2:binw:upperlim;
                            pdf1 = counts / sum(counts) / binw;
                            
                            plot(binCenters, pdf1, 'Color', color2(i,:), 'LineWidth', C.plot.lineWidth);
                            
                            % Add legend entry
                            text('Units', 'Normalized', 'Position', [legend_posX 0.75-i*0.12], ...
                                 'string', sprintf('%.2f', b1s(i)*0.1), 'Interpreter', 'latex', ...
                                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                                 'FontSize', C.fonts.annotationSize, 'color', color2(i,:));
                        else
                            warning('No valid weights found for b1 = %.2f', b1);
                        end
                    else
                        warning('Required fields not found in weight file for b1 = %.2f', b1);
                    end
                else
                    warning('Weight file not found for b1 = %.2f', b1);
                end
                
            catch loadError
                warning('Error loading weight data for b1 = %.2f: %s', b1, loadError.message);
            end
        end
        
        % Format plot
        set(gca, 'XScale', 'log');
        xlim([binw/2 upperlim]);
        set(gca, 'XTick', [1e-1, 1, 10, 100], 'XTickLabel', {'10^{-1}', '10^0', '10^1', '10^2'});
        
        % Add labels
        text('Units', 'Normalized', 'Position', [0.5 -0.2], ...
             'string', '$$\log(J_{\mathrm{scale}})$$', 'Interpreter', 'latex', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.labelSize, 'color', 'k');
        
        text('Units', 'Normalized', 'Position', [legend_posX 0.78], ...
             'string', '$$\mathrm{Var}(J_{\mathrm{scale}})$$', 'Interpreter', 'latex', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', 'k');
        
    catch ME
        warning('Error processing Panel A: %s', ME.message);
        plot(1:10, rand(1,10), 'k'); % Placeholder
    end
    
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Balance index and variance
    axes(AH(2));
    
    try
        mu1 = zeros(length(b1s)+1, 1);
        varI = zeros(length(b1s)+1, 1);
        
        % Load data for each b1 value
        for i = 1:length(b1s)
            b1 = b1s(i);
            filename = sprintf('Mean_Current_wide_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%s_bgamma_ex_%.1f_ix_%.1f_ee_%.1f_ie_%.1f_ei_%.1f_ii_%.1f', ...
                              1/scaleJ, 1/scaledown_ffwd, sigma_current, att_cond, sim_flag, ...
                              b1, b1, b1, b1, b1, b1);
            filename = strrep(filename, '.', 'd');
            
            try
                data1 = load(fullfile(rootFolder, 'data', [filename '.mat']));
                if isfield(data1, 'mu') && isfield(data1, 'I')
                    mu1(i) = data1.mu;
                    varI(i) = var(data1.I(1:C.network.Ne));
                else
                    warning('Required fields not found for b1 = %.2f', b1);
                    mu1(i) = NaN;
                    varI(i) = NaN;
                end
            catch
                warning('Could not load balance data for b1 = %.2f', b1);
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
            plot(1:5, mu1, 'Color', coloraxis(1,:), 'LineWidth', C.plot.lineWidth);
            for i = 1:5
                if validMu(i)
                    scatter(i, mu1(i), C.plot.scatterPointSize, color2(i,:), 'filled');
                end
            end
        end
        ylabel('balance index', 'FontSize', C.fonts.labelSize);
        ylim([0 0.7]);
        set(gca, 'YColor', coloraxis(1,:));
        
        yyaxis right;
        validVar = ~isnan(varI);
        if any(validVar)
            plot(1:5, varI, 'Color', coloraxis(2,:), 'LineWidth', C.plot.lineWidth);
            for i = 1:5
                if validVar(i)
                    scatter(i, varI(i), C.plot.scatterPointSize, color2(i,:), 'filled');
                end
            end
        end
        ylabel('var(I_{tot})', 'FontSize', C.fonts.labelSize);
        ylim([0 0.35]);
        set(gca, 'YColor', coloraxis(2,:));
        
        xlim([1 5]);
        set(gca, 'XTick', 1:5, 'XTickLabel', {'0.02', '0.05', '0.5', '0.83', 'default'});
        
    catch ME
        warning('Error processing Panel B: %s', ME.message);
        plot(1:4, rand(1,4), 'k');
    end
    
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C: Normalization index distributions
    axes(AH(3)); % Note: Panel C is actually the 6th axes
    
    try
        for i = 1:length(b1s)
            b1 = b1s(i);
            filename = sprintf('Norm_ind_wide_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%s_bgamma_ex_%.1f_ix_%.1f_ee_%.1f_ie_%.1f_ei_%.1f_ii_%.1f', ...
                              1/scaleJ, 1/scaledown_ffwd, sigma_current, att_cond, sim_flag, ...
                              b1, b1, b1, b1, b1, b1);
            filename = strrep(filename, '.', 'd');
            
            try
                data1 = load(fullfile(rootFolder, 'data', [filename '.mat']), 'pdf1');
                if isfield(data1, 'pdf1') && ~isempty(data1.pdf1)
                    xvals = 0.05:0.1:5;
                    if length(data1.pdf1) == length(xvals)
                        plot(xvals, data1.pdf1, 'Color', color2(i,:), 'LineWidth', C.plot.lineWidth);
                    end
                end
            catch
                warning('Could not load normalization distribution for b1 = %.2f', b1);
            end
        end
        
        % Add default case
        try
            defaultData = load(fullfile(rootFolder, 'data', 'Fig1_NormIndDistribution.mat'));
            if isfield(defaultData, 'norm1')
                norm1_filtered = defaultData.norm1(defaultData.norm1 < 5);
                counts = histcounts(norm1_filtered, 0:0.1:5);
                binCenters = 0.05:0.1:4.95;
                pdf1 = counts / length(norm1_filtered) / 0.1;
                plot(binCenters, pdf1, 'Color', color2(4,:), 'LineWidth', C.plot.lineWidth);
            end
        catch
            warning('Could not load default normalization distribution');
        end
        
        % Add legend
        text('Position', [3.5 0.8], 'string', 'Default', 'Interpreter', 'latex', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', color2(4,:));
        plot([2.5 3.2], [1 1]*0.8, 'Color', 'k', 'LineStyle', '-', 'LineWidth', C.plot.lineWidth);
        
    catch ME
        warning('Error processing Panel C: %s', ME.message);
    end
    
    xlabel('normalization index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    xlim([0 5]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panels D1-D3: Correlation heatmaps
    iAs = [4,5,6,7]; % Panels D1, D2, D3
    
    for panelIdx = 1:length(b1s)
        axes(AH(iAs(panelIdx)));
        
        try
            b1 = b1s(panelIdx);
            filename = sprintf('Norm_rsc_wide_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%s_bgamma_ex_%.1f_ix_%.1f_ee_%.1f_ie_%.1f_ei_%.1f_ii_%.1f', ...
                              1/scaleJ, 1/scaledown_ffwd, sigma_current, att_cond, sim_flag, ...
                              b1, b1, b1, b1, b1, b1);
            filename = strrep(filename, '.', 'd');
            
            data1 = load(fullfile(rootFolder, 'data', [filename '.mat']));
            
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
                
                % Panel title
                text('Units', 'Normalized', 'Position', [0.5 1.12], ...
                     'string', sprintf('%.2f', b1s(panelIdx)*Jei), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', C.fonts.titleSize, 'color', color2(panelIdx,:));
                
            else
                warning('Invalid data structure for panel D%d', panelIdx);
                imagesc(rand(20, 20));
            end
            
        catch loadError
            warning('Could not load correlation data for panel D%d: %s', panelIdx, loadError.message);
            imagesc(rand(20, 20)); % Placeholder
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
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_BroadWeight', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_BroadWeight:GeneralError', ...
          'Failed to generate broad weight supplementary figure: %s', ME.message);
end

end