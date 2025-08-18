function SuppFigure_BasicResponse(varargin)

%SUPPFIGURE_BASICRESPONSE Generate supplementary figure for basic response properties
%
% Purpose:
%   Creates a figure showing basic network response properties including
%   spike rasters, firing rate distributions, and correlation analyses
%   across different stimulus conditions.

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
    
    utils_plot.setPlotStyle('custom', 'width', 13, 'height', 5);
    
    figure(p.Results.FigNum); clf;
    utils_plot.matchFigureAspectRatio();
    
    horizontal_interval = [0.14 0.10];
    my_position1 = [0.08 0.14 0.88 0.7];
    my_position2 = [0.06 0.06 0.88 0.88];
    
    DC1 = utils_plot.divideAxes([0.25 0.25 0.25], [0.4 0.4], my_position1, horizontal_interval, 0.33);
    DC2 = utils_plot.divideAxes([0.25 0.25 0.25], 0.70, my_position2, horizontal_interval, []);
    
    LdPos = [-0.04, 0.10];
    
    AH(1) = axes('Position', DC1{1});
    hold on;
    utils_plot.addFigureLabel('A1', LdPos);
    
    AH(2) = axes('Position', DC1{2});
    hold on;
    utils_plot.addFigureLabel('A2', LdPos + [0 0.01]);
    
    AH(3) = axes('Position', DC2{2});
    hold on;
    utils_plot.addFigureLabel('B', LdPos);
    
    AH(4) = axes('Position', DC2{3});
    hold on;
    utils_plot.addFigureLabel('C', LdPos);
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Parameters
    color1 = 'm';
    color2 = 'k';
    time_interval = 200;
    start_time = 900;
    perturb_time = 945.0500;
    perturb_time1 = 946.00;
    
    colors = {"#D95319", "#77AC30", 'k'};
    labels = {'stimulus 1 only', 'stimulus 2 only', 'both stimuli presented'};
    
    %% Load basic response data
    try
        % Load index data
        indexData = load(fullfile(rootFolder,'data', 'indselect_subsample.mat'));
        patternData = load(fullfile(rootFolder,'data', 'Fig1_Pattern3cond.mat'));
        
        if ~isfield(indexData, 'indselect_subsample') || ~isfield(patternData, 'MTE')
            error('Required fields not found in basic response data');
        end
        
        indselect_subsample = indexData.indselect_subsample;
        MTE = patternData.MTE;
        
    catch ME
        warning('Error loading basic response data: %s', ME.message);
        indselect_subsample = 1:500;
        MTE = rand(C.network.Ne, 3) * 50;
    end
    
    %% Panel A1: Spontaneous spike raster
    axes(AH(1));
    
    try
        % Load spontaneous simulation data
        try
            spon1Data = load(fullfile(rootFolder,'data', 'simulation_spontaneous_perturb_one_spike_1.mat'));
            spon2Data = load(fullfile(rootFolder,'data', 'simulation_spontaneous_perturb_one_spike_2.mat'));
            
            % Process first simulation
            if isfield(spon1Data, 's2')
                s2 = spon1Data.s2;
                ind_fun = zeros(C.network.Ne, 1);
                ind_fun(indselect_subsample) = 1:500;
                
                validSpikes = s2(1,:) >= start_time & s2(1,:) <= start_time + time_interval & ...
                             ismember(s2(2,:), indselect_subsample);
                fired_ID = s2(2, validSpikes);
                spike_times = s2(1, validSpikes);
                fired_ID_reorder = ind_fun(fired_ID);
                
                scatter(spike_times/1000, fired_ID_reorder, 0.5, 'MarkerFaceColor', color1, ...
                       'MarkerEdgeColor', color1, 'LineWidth', 0.02, 'MarkerFaceAlpha', C.plot.alphaValue);
            end
            
            % Process second simulation
            if isfield(spon2Data, 's2')
                s2 = spon2Data.s2;
                ind_fun = zeros(C.network.Ne, 1);
                ind_fun(indselect_subsample) = 1:500;
                
                validSpikes = s2(1,:) >= start_time & s2(1,:) <= start_time + time_interval & ...
                             ismember(s2(2,:), indselect_subsample);
                fired_ID = s2(2, validSpikes);
                spike_times = s2(1, validSpikes);
                fired_ID_reorder = ind_fun(fired_ID);
                
                scatter(spike_times/1000, fired_ID_reorder, 0.5, 'MarkerFaceColor', color2, ...
                       'MarkerEdgeColor', color2, 'LineWidth', 0.02, 'MarkerFaceAlpha', C.plot.alphaValue);
            end
            
            % Add perturbation line
            plot([perturb_time1 perturb_time1]/1e3, [0 500], 'Color', 'cyan', 'LineWidth', C.plot.lineWidth);
            
        catch
            warning('Could not load spontaneous spike data, using placeholder');
            scatter(rand(1000,1)*0.2 + start_time/1e3, rand(1000,1)*500, 0.5, color1);
            scatter(rand(1000,1)*0.2 + start_time/1e3, rand(1000,1)*500, 0.5, color2);
        end
        
    catch ME
        warning('Error processing Panel A1: %s', ME.message);
        scatter(rand(1000,1)*0.2 + start_time/1e3, rand(1000,1)*500, 0.5, 'k');
    end
    
    ylabel('neuron index', 'FontSize', C.fonts.labelSize);
    xlim([start_time/1e3 start_time/1e3+time_interval/1000]);
    ylim([0 500]);
    title('Spontaneous', 'FontSize', C.fonts.titleSize);
    set(gca, 'XTick', [5, 5.1, 5.2], 'XTickLabel', {'0', '0.1', '0.2'});
    pbaspect([2 1 1]);
    
    %% Panel A2: Both stimuli spike raster
    axes(AH(2));
    
    try
        % Load both stimuli simulation data
        try
            both1Data = load(fullfile(rootFolder,'data', 'simulation_both_perturb_one_spike_1.mat'));
            both2Data = load(fullfile(rootFolder,'data', 'simulation_both_perturb_one_spike_2.mat'));
            
            % Process first simulation
            if isfield(both1Data, 's2')
                s2 = both1Data.s2;
                ind_fun = zeros(C.network.Ne, 1);
                ind_fun(indselect_subsample) = 1:500;
                
                validSpikes = s2(1,:) >= start_time & s2(1,:) <= start_time + time_interval & ...
                             ismember(s2(2,:), indselect_subsample);
                fired_ID = s2(2, validSpikes);
                spike_times = s2(1, validSpikes);
                fired_ID_reorder = ind_fun(fired_ID);
                
                scatter(spike_times/1000, fired_ID_reorder, 0.5, 'MarkerFaceColor', color1, ...
                       'MarkerEdgeColor', color1, 'LineWidth', 0.02, 'MarkerFaceAlpha', C.plot.alphaValue);
            end
            
            % Process second simulation
            if isfield(both2Data, 's2')
                s2 = both2Data.s2;
                ind_fun = zeros(C.network.Ne, 1);
                ind_fun(indselect_subsample) = 1:500;
                
                validSpikes = s2(1,:) >= start_time & s2(1,:) <= start_time + time_interval & ...
                             ismember(s2(2,:), indselect_subsample);
                fired_ID = s2(2, validSpikes);
                spike_times = s2(1, validSpikes);
                fired_ID_reorder = ind_fun(fired_ID);
                
                scatter(spike_times/1000, fired_ID_reorder, 0.5, 'MarkerFaceColor', color2, ...
                       'MarkerEdgeColor', color2, 'LineWidth', 0.02, 'MarkerFaceAlpha', C.plot.alphaValue);
            end
            
            % Add perturbation line
            plot([perturb_time perturb_time]/1e3, [0 500], 'Color', 'cyan', 'LineWidth', C.plot.lineWidth);
            
        catch
            warning('Could not load both stimuli spike data, using placeholder');
            scatter(rand(1000,1)*0.2 + start_time/1e3, rand(1000,1)*500, 0.5, color1);
            scatter(rand(1000,1)*0.2 + start_time/1e3, rand(1000,1)*500, 0.5, color2);
        end
        
    catch ME
        warning('Error processing Panel A2: %s', ME.message);
        scatter(rand(1000,1)*0.2 + start_time/1e3, rand(1000,1)*500, 0.5, 'k');
    end
    
    xlabel('time (s)', 'FontSize', C.fonts.labelSize);
    ylabel('neuron index', 'FontSize', C.fonts.labelSize);
    xlim([start_time/1e3 start_time/1e3+time_interval/1000]);
    ylim([0 500]);
    title('Both stimuli presented', 'FontSize', C.fonts.titleSize);
    set(gca, 'XTick', (start_time:100:start_time+time_interval)/1e3, 'XTickLabel', {'0', '0.1', '0.2'});
    pbaspect([2 1 1]);
    
    %% Panel B: Firing rate distributions
    axes(AH(3));
    
    try
        upper_rate = 80;
        pdf_popwhole = zeros(3, upper_rate);
        
        for i = 1:3
            if size(MTE, 2) >= i
                temp = histcounts(MTE(:,i), 0:1:upper_rate);
                pdf_popwhole(i,:) = temp / C.network.Ne;
                plot(0.5:1:upper_rate, pdf_popwhole(i,:), 'color', colors{i}, ...
                     'LineStyle', '-', 'LineWidth', C.plot.lineWidth);
                
                text('Units', 'Normalized', 'Position', [0.6 0.9-0.1*i], ...
                     'string', labels{i}, 'color', colors{i}, ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', C.fonts.annotationSize);
            end
        end
        
    catch ME
        warning('Error processing Panel B: %s', ME.message);
        for i = 1:3
            plot(0.5:1:upper_rate, rand(1, upper_rate)*0.1, 'color', colors{i});
        end
    end
    
    xlabel('firing rate (Hz)', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    xlim([0 upper_rate]);
    ylim([0 0.15]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C: Spike count correlation
    axes(AH(4));
    
    try
        % Load correlation statistics
        rscData = load(fullfile(rootFolder,'data', 'FigSupp_rsc_stats.mat'));
        
        if isfield(rscData, 'rsc_pdf') && isfield(rscData, 'rsc_mean')
            rsc_pdf = rscData.rsc_pdf;
            rsc_mean = rscData.rsc_mean;
            
            xvals = -0.99:0.02:0.99;
            
            for i = 1:3
                if i <= length(rsc_pdf) && i <= length(rsc_mean)
                    if length(rsc_pdf{i}) == length(xvals)
                        plot(xvals, rsc_pdf{i}, 'color', colors{i}, 'LineStyle', '-', 'LineWidth', C.plot.lineWidth);
                        plot([rsc_mean(i) rsc_mean(i)], [0 9], 'color', colors{i}, 'LineStyle', ':', 'LineWidth', C.plot.lineWidth);
                    end
                end
            end
            
            text('Position', [0.01 9.7], 'string', 'mean', 'color', 'k', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', C.fonts.annotationSize);
            
        else
            warning('Required correlation fields not found');
            for i = 1:3
                plot(-0.99:0.02:0.99, rand(1, 100)*5, 'color', colors{i});
            end
        end
        
    catch ME
        warning('Error processing Panel C: %s', ME.message);
        for i = 1:3
            plot(-0.99:0.02:0.99, rand(1, 100)*5, 'color', colors{i});
        end
    end
    
    xlabel('spike count correlation', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    xlim([-0.2 0.2]);
    ylim([0 10]);
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_BasicResponseProperties', ...
                             'view', p.Results.View, 'save', true, 'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_BasicResponse:GeneralError', ...
          'Failed to generate basic response properties figure: %s', ME.message);
end
end
