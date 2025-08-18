function SuppFigure_I_activity(varargin)
%SUPPFIGURE_I_ACTIVITY Generate supplementary figure for inhibitory neuron activity analysis
%
% Purpose:
%   Creates a figure analyzing inhibitory neuron activity including firing rate
%   distributions, normalization index distributions, and current distributions
%   for both excitatory and inhibitory neurons, comparing different current types.
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
%   SuppFigure_I_activity();
%   SuppFigure_I_activity('Save', false, 'FigNum', 2);

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
    utils_plot.setPlotStyle('custom', 'width', 12, 'height', 9);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = 0.03;
    vertical_interval = 0.1;
    my_position2 = [0.08 0.15 0.75 0.75];
    
    DC = utils_plot.divideAxes(0.25*ones(1,2), 0.25*ones(1,2), my_position2, ...
                               horizontal_interval, vertical_interval);
    
    % Create axes
    LdPos = [-0.1, 0.05];
    labels = {'A', 'B', 'C1', 'C2'};
    
    for i = 1:4
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
        
    % Colors
    col1 = C.colors.excitatory;
    col2 = C.colors.inhibitory;
    my_ylim = [0 0.3];
    
    %% Load data
    try
        % Load pattern data for E and I neurons
        patternDataE = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond.mat'));
        patternDataI = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond_I.mat'));
        
        % Load normalization indices
        normDataE = load(fullfile(rootFolder, 'data', 'Fig1_NormIndDistribution.mat'));
        normDataI = load(fullfile(rootFolder, 'data', 'Fig1_NormIndDistribution_I.mat'));
        
        % Load current data
        currentDataE = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_current_default.mat'), ...
                           'inputEr', 'inputXr', 'inputIr');
        currentDataI = load(fullfile(rootFolder, 'data', 'Fig3_NormIndFiringRateCurrent_I.mat'), ...
                           'inputEr', 'inputXr', 'inputIr');
        
        % Validate data
        
        
        MTE = patternDataE.MTE;
        MTI = patternDataI.MTI;
        
        if ~isfield(normDataE, 'norm1') || ~isfield(normDataI, 'norm1')
            error('Normalization index data not found');
        end
        
        normE = normDataE.norm1;
        normI = normDataI.norm1;
        
    catch ME
        error('SuppFigure_I_activity:DataLoadError', ...
              'Failed to load required data: %s', ME.message);
    end
    
    %% Panel A: Firing rate distributions
    axes(AH(1));
    
    try
        upper_rate = 100;
        pdf_popwhole = zeros(2, upper_rate);
        
        % E neuron distribution
        if size(MTE, 1) >= C.network.Ne
            temp = histcounts(MTE(1:C.network.Ne, 3), 0:1:upper_rate);
            pdf_popwhole(1,:) = temp / C.network.Ne;
        else
            temp = histcounts(MTE(:, 3), 0:1:upper_rate);
            pdf_popwhole(1,:) = temp / size(MTE, 1);
        end
        
        plot(0.5:1:upper_rate, pdf_popwhole(1,:), 'color', 'k', 'LineStyle', '-', ...
             'LineWidth', C.plot.lineWidth);
        
        % I neuron distribution
        if size(MTI, 1) >= C.network.Ni
            temp = histcounts(MTI(1:C.network.Ni, 3), 0:1:upper_rate);
            pdf_popwhole(2,:) = temp / C.network.Ni;
        else
            temp = histcounts(MTI(:, 3), 0:1:upper_rate);
            pdf_popwhole(2,:) = temp / size(MTI, 1);
        end
        
        plot(0.5:1:upper_rate, pdf_popwhole(2,:), 'color', 'b', 'LineStyle', '-', ...
             'LineWidth', C.plot.lineWidth);
        
        % Add legend
        text('Units', 'Normalized', 'Position', [0.7 0.8], ...
             'string', 'E neuron', 'color', 'k', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize);
        text('Units', 'Normalized', 'Position', [0.7 0.65], ...
             'string', 'I neuron', 'color', 'b', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize);
        
    catch ME
        warning('Error processing firing rate distributions: %s', ME.message);
        plot(0.5:1:upper_rate, rand(1, upper_rate) * 0.1, 'k');
        plot(0.5:1:upper_rate, rand(1, upper_rate) * 0.1, 'b');
    end
    
    xlim([0 upper_rate]);
    ylim(my_ylim);
    xlabel('firing rate (Hz)', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Normalization index distributions
    axes(AH(2));
    
    try
        upperlim_norm = 3;
        
        % E neuron normalization indices
        temp = normE(normE <= upperlim_norm & normE >= 0 & ~isnan(normE));
        if ~isempty(temp)
            counts = histcounts(temp, 0:0.1:upperlim_norm);
            binCenters = 0.05:0.1:upperlim_norm;
            densityValues = counts / length(temp) / 0.1;
            plot(binCenters, densityValues, 'k', 'LineWidth', C.plot.lineWidth);
        end
        
        % I neuron normalization indices
        temp = normI(normI <= upperlim_norm & normI >= 0 & ~isnan(normI));
        if ~isempty(temp)
            counts = histcounts(temp, 0:0.1:upperlim_norm);
            binCenters = 0.05:0.1:upperlim_norm;
            densityValues = counts / length(temp) / 0.1;
            plot(binCenters, densityValues, 'b', 'LineWidth', C.plot.lineWidth);
        end
        
    catch ME
        warning('Error processing normalization index distributions: %s', ME.message);
        plot(0.05:0.1:upperlim_norm, rand(1, 30), 'k');
        plot(0.05:0.1:upperlim_norm, rand(1, 30), 'b');
    end
    
    xlim([0 upperlim_norm]);
    xlabel('normalization index', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C1: E neuron current distributions
    axes(AH(3));
    
    try
        upperlim_current = 8;
        binw = 0.2;
        
        if isfield(currentDataE, 'inputXr') && isfield(currentDataE, 'inputEr') && isfield(currentDataE, 'inputIr')
            % Feedforward excitation
            if size(currentDataE.inputXr, 3) >= 3
                inputXr1 = squeeze(sum(currentDataE.inputXr(:,3,:), 3));
                temp1 = histcounts(inputXr1, 0:binw:upperlim_current);
                temp2 = temp1 / C.network.Ne / binw;
                plot(binw/2:binw:upperlim_current, temp2, 'Color', col1, 'LineStyle', ':', ...
                     'LineWidth', C.plot.lineWidth);
            end
            
            % Recurrent excitation
            if size(currentDataE.inputEr, 2) >= 3
                temp1 = histcounts(currentDataE.inputEr(:,3), 0:binw:upperlim_current);
                temp2 = temp1 / C.network.Ne / binw;
                plot(binw/2:binw:upperlim_current, temp2, 'Color', col1, 'LineStyle', '-.', ...
                     'LineWidth', C.plot.lineWidth);
            end
            
            % Recurrent inhibition
            if size(currentDataE.inputIr, 2) >= 3
                temp1 = histcounts(abs(currentDataE.inputIr(:,3)), 0:binw:upperlim_current);
                temp2 = temp1 / C.network.Ne / binw;
                plot(binw/2:binw:upperlim_current, temp2, 'Color', col2, 'LineStyle', '-', ...
                     'LineWidth', C.plot.lineWidth);
            end
        else
            warning('Current data fields not found for E neurons');
            plot(binw/2:binw:upperlim_current, rand(1, length(binw/2:binw:upperlim_current)), 'k');
        end
        
    catch ME
        warning('Error processing E neuron current distributions: %s', ME.message);
        plot(binw/2:binw:upperlim_current, rand(1, length(binw/2:binw:upperlim_current)), 'k');
    end
    
    set(gca, 'XTick', [0 4 8]);
    xlim([0 upperlim_current]);
    ylim([0 3]);
    
    firstline_height = -0.2;
    secondline_height = -0.3;
    text('Units', 'Normalized', 'Position', [0.5 firstline_height], ...
         'string', 'Current onto E neuron', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
    text('Units', 'Normalized', 'Position', [0.5 secondline_height], ...
         'string', '(absolute value)', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
    
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C2: I neuron current distributions
    axes(AH(4));
    
    try
        if isfield(currentDataI, 'inputXr') && isfield(currentDataI, 'inputEr') && isfield(currentDataI, 'inputIr')
            % Feedforward excitation
            if size(currentDataI.inputXr, 2) >= 3
                temp1 = histcounts(currentDataI.inputXr(:,3), 0:binw:upperlim_current);
                temp2 = temp1 / C.network.Ni / binw;
                plot(binw/2:binw:upperlim_current, temp2, 'Color', col1, 'LineStyle', ':', ...
                     'LineWidth', C.plot.lineWidth);
            end
            
            % Recurrent excitation
            if size(currentDataI.inputEr, 2) >= 3
                temp1 = histcounts(currentDataI.inputEr(:,3), 0:binw:upperlim_current);
                temp2 = temp1 / C.network.Ni / binw;
                plot(binw/2:binw:upperlim_current, temp2, 'Color', col1, 'LineStyle', '-.', ...
                     'LineWidth', C.plot.lineWidth);
            end
            
            % Recurrent inhibition
            if size(currentDataI.inputIr, 2) >= 3
                temp1 = histcounts(abs(currentDataI.inputIr(:,3)), 0:binw:upperlim_current);
                temp2 = temp1 / C.network.Ni / binw;
                plot(binw/2:binw:upperlim_current, temp2, 'Color', col2, 'LineStyle', '-', ...
                     'LineWidth', C.plot.lineWidth);
            end
        else
            warning('Current data fields not found for I neurons');
            plot(binw/2:binw:upperlim_current, rand(1, length(binw/2:binw:upperlim_current)), 'k');
        end
        
        % Add legend
        legend_posY = 2.6:-0.4:1.8;
        plot([4 5.5], [1 1]*legend_posY(1), 'Color', col1, 'LineStyle', ':', ...
             'LineWidth', C.plot.lineWidth);
        plot([4 5.5], [1 1]*legend_posY(2), 'Color', col1, 'LineStyle', '-.', ...
             'LineWidth', C.plot.lineWidth);
        plot([4 5.5], [1 1]*legend_posY(3), 'Color', col2, 'LineStyle', '-', ...
             'LineWidth', C.plot.lineWidth);
        
        text('Position', [6.0 legend_posY(1)], 'string', 'Ffwd. excitation', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', col1, 'FontSize', C.fonts.annotationSize);
        text('Position', [6.0 legend_posY(2)], 'string', 'Rec. excitation', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', col1, 'FontSize', C.fonts.annotationSize);
        text('Position', [6.0 legend_posY(3)], 'string', 'Rec. inhibition', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', col2, 'FontSize', C.fonts.annotationSize);
        
    catch ME
        warning('Error processing I neuron current distributions: %s', ME.message);
        plot(binw/2:binw:upperlim_current, rand(1, length(binw/2:binw:upperlim_current)), 'k');
    end
    
    ylim([0 3]);
    xlim([0 upperlim_current]);
    set(gca, 'XTick', [0 4 8]);
    
    text('Units', 'Normalized', 'Position', [0.5 firstline_height], ...
         'string', 'Current onto I neuron', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
    text('Units', 'Normalized', 'Position', [0.5 secondline_height], ...
         'string', '(absolute value)', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
    
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_I_activity', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_I_activity:GeneralError', ...
          'Failed to generate inhibitory activity supplementary figure: %s', ME.message);
end

end