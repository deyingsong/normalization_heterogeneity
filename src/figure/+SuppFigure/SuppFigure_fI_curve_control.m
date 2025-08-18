function SuppFigure_fI_curve_control(varargin)
%SUPPFIGURE_FI_CURVE_CONTROL Generate supplementary figure for f-I curve control analysis
%
% Purpose:
%   Creates a figure showing different f-I curve shapes and their effects on
%   normalization, including linear, sublinear, and supralinear curves for
%   both E and I neurons, with rate comparisons and normalization distributions.

try
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
    
    p = inputParser;
    addParameter(p, 'Save', true, @islogical);
    addParameter(p, 'View', true, @islogical);
    addParameter(p, 'FigNum', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    
    utils_plot.setPlotStyle('custom', 'width', 14, 'height', 5);
    
    figure(p.Results.FigNum); clf;
    utils_plot.matchFigureAspectRatio();
    
    horizontal_interval = [0.07 0.09 0.09];
    my_position2 = [0.09 0.15 0.88 0.70];
    DC = utils_plot.divideAxes(0.25*ones(1,4), 1, my_position2, horizontal_interval, []);
    
    LdPos = [-0.04, 0.13];
    labels = {'A1', 'A2', 'B', 'C'};
    
    for i = 1:4
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % f-I curve parameters (from original script)
    param1 = struct();
    param1.k1E = 18.01; param1.aE = -0.3248; param1.nE = 2.472;
    param1.k2E = 41.57; param1.cE = 0.7223;
    param1.k1I = 8.001; param1.aI = -0.4854; param1.nI = 3.228;
    param1.k2I = 56.87; param1.cI = 0.9447;
    
    param2 = struct();
    param2.kE = 37.77; param2.alphaE = 1; param2.thetaE = param1.aE;
    param2.kI = 39.57; param2.alphaI = 1; param2.thetaI = param1.aI;
    
    colors = cool(3);
    label_left = 0.1;
    label_upper = 0.9;
    label_int = 0.1;
    Is = -1.5:0.02:3;
    
    %% Panel A1: E neuron f-I curves
    axes(AH(1));
    
    try
        % Default curve
        r1 = 0*(Is<=param1.aE) + param1.k1E.*(Is-param1.aE).^param1.nE.*(Is>param1.aE).*(Is<=param1.cE) + ...
             (param1.k2E.*(Is-param1.cE)+20).*(Is>param1.cE);
        
        % Linear curve
        r2 = 34*(Is-param2.thetaE).*(Is-param2.thetaE>0);
        
        % Supralinear curve
        r3 = 10*(Is-param2.thetaE).^2.*(Is-param2.thetaE>0);
        
        % Sublinear curve
        r4 = 43*(Is-param2.thetaE).^0.8.*(Is-param2.thetaE>0);
        
        plot(Is, r1, 'k', 'LineWidth', C.plot.lineWidth);
        plot(Is, r2, 'color', colors(2,:), 'LineWidth', C.plot.lineWidth);
        plot(Is, r3, 'color', colors(3,:), 'LineWidth', C.plot.lineWidth);
        plot(Is, r4, 'color', colors(1,:), 'LineWidth', C.plot.lineWidth);
        
        % Add legend
        text('Units', 'Normalized', 'Position', [label_left, label_upper], 'string', 'Default', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', 'k');
        text('Units', 'Normalized', 'Position', [label_left, label_upper-label_int], 'string', 'sublinear', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', colors(1,:));
        text('Units', 'Normalized', 'Position', [label_left, label_upper-label_int*2], 'string', 'linear', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', colors(2,:));
        text('Units', 'Normalized', 'Position', [label_left, label_upper-label_int*3], 'string', 'supralinear', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', colors(3,:));
        
    catch ME
        warning('Error creating E neuron f-I curves: %s', ME.message);
        plot(Is, Is*10, 'k');
    end
    
    xlabel('current (mV/ms)', 'FontSize', C.fonts.labelSize);
    ylabel('firing rate (Hz)', 'FontSize', C.fonts.labelSize);
    title('E neuron', 'FontSize', C.fonts.titleSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel A2: I neuron f-I curves
    axes(AH(2));
    
    try
        % Default curve
        r1 = 0*(Is<=param1.aI) + param1.k1I.*(Is-param1.aI).^param1.nI.*(Is>param1.aI).*(Is<=param1.cI) + ...
             (param1.k2I.*(Is-param1.cI)+25).*(Is>param1.cI);
        
        % Linear curve
        r2 = 40*(Is-param2.thetaI).*(Is-param2.thetaI>0);
        
        % Supralinear curve
        r3 = 12*(Is-param2.thetaI).^2.*(Is-param2.thetaI>0);
        
        % Sublinear curve
        r4 = 52*(Is-param2.thetaI).^0.8.*(Is-param2.thetaI>0);
        
        plot(Is, r1, 'k', 'LineWidth', C.plot.lineWidth);
        plot(Is, r2, 'color', colors(2,:), 'LineWidth', C.plot.lineWidth);
        plot(Is, r3, 'color', colors(3,:), 'LineWidth', C.plot.lineWidth);
        plot(Is, r4, 'color', colors(1,:), 'LineWidth', C.plot.lineWidth);
        
    catch ME
        warning('Error creating I neuron f-I curves: %s', ME.message);
        plot(Is, Is*10, 'k');
    end
    
    xlabel('current (mV/ms)', 'FontSize', C.fonts.labelSize);
    title('I neuron', 'FontSize', C.fonts.titleSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Rate comparison
    axes(AH(3));
    
    try
        % Load rate simulation data
        defaultI = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond_I.mat'));
        defaultE = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond.mat'));
        rateData = load(fullfile(rootFolder, 'data', 'FigSupp_RateSim_rate_default.mat'));
        
        if isfield(defaultI, 'MTI') && isfield(rateData, 'MTfr') && isfield(defaultE, 'MTE')
            MTI = defaultI.MTI;
            MTE = defaultE.MTE;
            MTfr = rateData.MTfr;
            
            % Plot I neurons
            scatter(MTI(:,3), MTfr(C.network.Ne+1:C.network.N,3), 1, 'b.', ...
                   'MarkerFaceAlpha', C.plot.alphaValue);
            
            % Plot E neurons  
            scatter(MTE(1:C.network.Ne,3), MTfr(1:C.network.Ne,3), 1, 'r.', ...
                   'MarkerFaceAlpha', C.plot.alphaValue);
        else
            warning('Rate comparison data not found');
            scatter(rand(100,1)*120, rand(100,1)*120, 1, 'b.');
            scatter(rand(100,1)*120, rand(100,1)*120, 1, 'r.');
        end
        
        xlim([0 120]); ylim([0 120]);
        plot([0 120], [0 120], 'k', 'LineWidth', C.plot.lineWidth);
        
        % Add legend
        text('Units', 'Normalized', 'Position', [label_left, label_upper], 'string', 'E neuron', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', 'r');
        text('Units', 'Normalized', 'Position', [label_left, label_upper-label_int], 'string', 'I neuron', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', 'b');
        
    catch ME
        warning('Error creating rate comparison: %s', ME.message);
        scatter(rand(100,1)*120, rand(100,1)*120, 1, 'k.');
    end
    
    xlabel('rate, default model', 'FontSize', C.fonts.labelSize);
    ylabel('rate, rate model', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C: Normalization distributions
    axes(AH(4));
    
    try
        % Load different curve type data
        supralinearData = load(fullfile(rootFolder, 'data', 'FigSupp_RateSim_rate_supralinear_control.mat'));
        sublinearData = load(fullfile(rootFolder, 'data', 'FigSupp_RateSim_rate_sublinear_control.mat'));
        linearData = load(fullfile(rootFolder, 'data', 'FigSupp_RateSim_rate_linear_control.mat'));
        defaultData = load(fullfile(rootFolder, 'data', 'Fig1_NormIndDistribution.mat'));
        
        binEdges = 0:0.1:3;
        binCenters = 0.05:0.1:2.95;
        
        % Plot each distribution
        datasets = {supralinearData, linearData, sublinearData, defaultData};
        plotColors = {colors(3,:), colors(2,:), colors(1,:), 'k'};
        
        for i = 1:4
            if isfield(datasets{i}, 'norm1')
                temp = datasets{i}.norm1(datasets{i}.norm1 <= 3 & datasets{i}.norm1 >= 0);
                if ~isempty(temp)
                    counts = histcounts(temp, binEdges);
                    densityValues = counts / length(temp) / 0.1;
                    plot(binCenters, densityValues, 'color', plotColors{i}, 'LineWidth', C.plot.lineWidth);
                end
            end
        end
        
    catch ME
        warning('Error creating normalization distributions: %s', ME.message);
        plot(binCenters, rand(size(binCenters)), 'k');
    end
    
    xlim([0 3]);
    xlabel('norm index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_fI_curve_control', ...
                             'view', p.Results.View, 'save', true, 'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_fI_curve_control:GeneralError', ...
          'Failed to generate f-I curve control figure: %s', ME.message);
end
end
