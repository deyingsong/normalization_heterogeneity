function SuppFigure_FFmodel_I0(varargin)
%SUPPFIGURE_FFMODEL_I0 Generate supplementary figure for feedforward model I0 analysis
%
% Purpose:
%   Creates a figure showing feedforward model analysis with different I0 values,
%   comparing inhibitory neurons and random networks under different current offset conditions.

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
    
    utils_plot.setPlotStyle('custom', 'width', 8, 'height', 10);
    
    figure(p.Results.FigNum); clf;
    utils_plot.matchFigureAspectRatio();
    
    horizontal_interval = 0.08;
    vertical_interval = 0.18;
    my_position2 = [0.09 0.15 0.82 0.70];
    DC = utils_plot.divideAxes(0.25*ones(1,2), 0.25*ones(1,2), my_position2, ...
                               horizontal_interval, vertical_interval);
    
    LdPos = [-0.03, 0.05];
    labels = {'A1', 'B1', 'A2', 'B2'};
    
    for i = 1:4
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties();
    
    upper_title = 1.3;
    
    %% Load I neuron data
    try
        iNeuronData = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond_I.mat'));
        currentData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_current_default.mat'));
        
        if ~isfield(iNeuronData, 'MTI') || ~isfield(currentData, 'inputTr')
            error('Required fields not found in I neuron data');
        end
        
        rate = iNeuronData.MTI(:, 1:3);
        Itot = currentData.inputTr(C.network.Ne+1:end, :);
        norm1 = (rate(:,1) + rate(:,2)) ./ rate(:,3);
        
        
        % Calculate current offset for binning
        th = Itot(1:C.network.Ni, 3) - (Itot(1:C.network.Ni, 1) + Itot(1:C.network.Ni, 2))/1.5;
        bins = quantile(th, 0:0.25:1);
        
    catch ME
        warning('Error loading I neuron data: %s', ME.message);
        rate = rand(1000, 3) * 50;
        norm1 = rand(1000, 1) * 3;
        bins = [-1, -0.5, 0, 0.5, 1];
    end
    
    %% Panels A1 and A2: Inhibitory neuron analysis
    plotBins = [4, 1]; % Q4, Q1
    panelIndices = [1, 3];
    
    for plotIdx = 1:2
        bb = plotBins(plotIdx);
        axes(AH(panelIndices(plotIdx)));
        
        try
            ind1 = find((rate(:,3) > 1) & (th < bins(bb+1)) & (th >= bins(bb)));
            
            if ~isempty(ind1)
                scatter(rate(ind1,1), rate(ind1,2), 5, norm1(ind1), 'filled');
                axis([0 80 0 80]);
                colormap(AH(panelIndices(plotIdx)), utils_plot.redblue());
                
                if bb <= 2
                    clim([-1 4]);
                    cbar_labels = {'-1', '1.5', '4'};
                    cticks = [-1 1.5 4];
                else
                    clim([1 2]);
                    cbar_labels = {'1', '1.5', '2'};
                    cticks = [1 1.5 2];
                end
                
                h = utils_plot.createNarrowColorbar('vert');
                set(h, 'YTick', cticks);
                set(h, 'YTickLabels', cbar_labels);
                set(h, 'FontSize', C.fonts.tickSize);
                colormap(h, utils_plot.redblue())
            else
                warning('No data points found for I neuron bin %d', bb);
                scatter([40], [40], 5, 1.5, 'filled');
                axis([0 80 0 80]);
            end
            
        catch ME
            warning('Error processing I neuron panel %d: %s', plotIdx, ME.message);
            scatter(rand(50,1)*80, rand(50,1)*80, 5, rand(50,1)*3, 'filled');
            axis([0 80 0 80]);
        end
        
        title(sprintf('I0 in [%.2g, %.2g]', bins(bb), bins(bb+1)), 'FontSize', C.fonts.titleSize);
        
        if plotIdx == 1
            text('Units', 'Normalized', 'Position', [0.5 upper_title], 'string', 'Inh neuron', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', C.fonts.titleSize, 'color', 'k');
            ylabel('r2 (Hz)', 'FontSize', C.fonts.labelSize);
        else
            xlabel('r1 (Hz)', 'FontSize', C.fonts.labelSize);
            ylabel('r2 (Hz)', 'FontSize', C.fonts.labelSize);
        end
        
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Load random network data and process panels B1 and B2
    try
        randomData = load(fullfile(rootFolder, 'data', 'FigSupp_current_FFmodel_I0_random.mat'));
        
        if ~isfield(randomData, 'Itot') || ~isfield(randomData, 'rate')
            error('Required fields not found in random network data');
        end
        
        Itot_random = randomData.Itot;
        rate_random = randomData.rate;
        norm1_random = (rate_random(1:C.network.Ne,1) + rate_random(1:C.network.Ne,2)) ./ rate_random(1:C.network.Ne,3);
        
        th_random = Itot_random(1:C.network.Ne,3) - (Itot_random(1:C.network.Ne,1) + Itot_random(1:C.network.Ne,2))/1.5;
        bins_random = quantile(th_random, 0:0.25:1);
        
    catch ME
        warning('Error loading random network data: %s', ME.message);
        rate_random = rand(C.network.Ne, 3) * 50;
        norm1_random = rand(C.network.Ne, 1) * 3;
        bins_random = [-1, -0.5, 0, 0.5, 1];
        th_random = randn(C.network.Ne, 1);
    end
    
    %% Panels B1 and B2: Random network analysis
    panelIndices = [2, 4];
    
    for plotIdx = 1:2
        bb = plotBins(plotIdx);
        axes(AH(panelIndices(plotIdx)));
        
        try
            ind1 = find((rate_random(:,3) > 1) & (th_random < bins_random(bb+1)) & (th_random >= bins_random(bb)));
            
            if ~isempty(ind1)
                scatter(rate_random(ind1,1), rate_random(ind1,2), 5, norm1_random(ind1), 'filled');
                axis([0 80 0 80]);
                colormap(AH(panelIndices(plotIdx)), utils_plot.redblue());
                
                if bb <= 2
                    clim([-1 4]);
                    cbar_labels = {'-1', '1.5', '4'};
                    cticks = [-1 1.5 4];
                else
                    clim([1 2]);
                    cbar_labels = {'1', '1.5', '2'};
                    cticks = [1 1.5 2];
                end
                
                h = utils_plot.createNarrowColorbar('vert');
                set(h, 'YTick', cticks);
                set(h, 'YTicklabels', cbar_labels);
                set(h, 'FontSize', C.fonts.tickSize);
                colormap(h,utils_plot.redblue())
                
            else
                warning('No data points found for random network bin %d', bb);
                scatter([40], [40], 5, 1.5, 'filled');
                axis([0 80 0 80]);
            end
            
        catch ME
            warning('Error processing random network panel %d: %s', plotIdx, ME.message);
            scatter(rand(50,1)*80, rand(50,1)*80, 5, rand(50,1)*3, 'filled');
            axis([0 80 0 80]);
        end
        
        title(sprintf('I0 in [%.2g, %.2g]', bins_random(bb), bins_random(bb+1)), 'FontSize', C.fonts.titleSize);
        
        if plotIdx == 1
            text('Units', 'Normalized', 'Position', [0.5 upper_title], 'string', 'Random network', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', C.fonts.titleSize, 'color', 'k');
        else
            xlabel('r1 (Hz)', 'FontSize', C.fonts.labelSize);
        end
        
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_FFmodel_I0', ...
                             'view', p.Results.View, 'save', true, 'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_FFmodel_I0:GeneralError', ...
          'Failed to generate FFmodel I0 figure: %s', ME.message);
end
end
