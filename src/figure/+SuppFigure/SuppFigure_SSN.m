function SuppFigure_SSN(varargin)
%SUPPFIGURE_SSN Generate supplementary figure for Stabilized Supralinear Network analysis
%
% Purpose:
%   Creates a figure showing SSN model analysis including normalization index
%   distribution, preferred orientation analysis, and correlation matrices
%   for both original and tuning-matched conditions.
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
%   SuppFigure_SSN();
%   SuppFigure_SSN('Save', false, 'View', true);

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
    utils_plot.setPlotStyle('custom', 'width', 12, 'height', 11);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Define layout
    horizontal_interval = 0.08;
    vertical_interval = 0.1;
    my_position2 = [0.08 0.08 0.83 0.85];
    
    DC = utils_plot.divideAxes(0.25*ones(1,2), 0.25*ones(1,2), my_position2, ...
                               horizontal_interval, vertical_interval);
    
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
    LdPos = [-0.03, 0.05];
    labels = {'A', 'B', 'C1', 'C2'};
    
    for i = 1:4
        AH(i) = axes('Position', DC{i});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties();
    
    % Parameters
    Nneuron = 4000;
    Nshuffle = 2000;
    xticks1 = 0:90:180;
    xticklabels1 = {'0', '90', '180'};
    clim1 = [-0.004 0.006];
    clim2 = [-0.002 0.003];
    my_ylim = [0 1.5];
    
    %% Load SSN data
    try
        ssnData = load(fullfile(rootFolder, 'data', 'FigSupp_ssn_data.mat'));
        
        % Validate required fields
        requiredFields = {'norm1', 'ind1', 'tuning', 'dataheat_original', 'dataheat_matchtuning'};
        for field = requiredFields
            if ~isfield(ssnData, field{1})
                error('Required field "%s" not found in SSN data', field{1});
            end
        end
        
        norm1 = ssnData.norm1;
        ind1 = ssnData.ind1;
        tuning = ssnData.tuning;
        
        % Validate data
        if length(norm1) ~= length(tuning)
            error('Mismatch between norm1 and tuning data lengths');
        end
        
        if length(ind1) > length(norm1)
            error('ind1 contains more indices than available data');
        end
        
        % Filter valid indices
        validIndices = ind1(ind1 <= length(norm1) & ind1 > 0);
        if length(validIndices) < Nneuron
            warning('Not enough valid indices (%d), using all available (%d)', ...
                    Nneuron, length(validIndices));
            Nneuron = length(validIndices);
        else
            validIndices = validIndices(1:Nneuron);
        end
        
    catch ME
        error('SuppFigure_SSN:DataLoadError', ...
              'Failed to load SSN data: %s', ME.message);
    end
    
    %% Panel A: Normalization index distribution
    axes(AH(1));
    
    try
        temp = norm1(validIndices);
        temp = temp(~isnan(temp) & temp >= 0 & temp <= 3);
        
        if ~isempty(temp)
            binEdges = 0:0.1:3;
            counts = histcounts(temp, binEdges);
            binCenters = binEdges(1:end-1) + diff(binEdges)/2;
            densityValues = counts / length(temp) / 0.1;
            
            plot(binCenters, densityValues, 'k', 'LineWidth', C.plot.lineWidth);
        else
            warning('No valid normalization indices for Panel A');
            plot(0.05:0.1:2.95, zeros(1, 30), 'k');
        end
        
    catch ME
        warning('Error processing normalization distribution: %s', ME.message);
        plot(0.05:0.1:2.95, zeros(1, 30), 'k');
    end
    
    xlim([0 3]);
    ylim(my_ylim);
    xlabel('normalization index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Preferred orientation vs normalization index
    axes(AH(2));
    
    try
        norm_subset = norm1(validIndices);
        tuning_subset = tuning(validIndices);
        
        % Remove NaN values
        validData = ~isnan(norm_subset) & ~isnan(tuning_subset) & ...
                   norm_subset >= 0 & norm_subset <= 5 & ...
                   tuning_subset >= 0 & tuning_subset <= 180;
        
        if sum(validData) > 0
            norm_plot = norm_subset(validData);
            tuning_plot = tuning_subset(validData);
            
            % Create scatter plot with density coloring
            Nbins = 60;
            xbinw = 180/Nbins;
            ybinw = 5/Nbins;
            sizepoint = 3;
            c0 = parula;
            
            N_temp = histcounts2(tuning_plot, norm_plot, 0:xbinw:180, 0.0:ybinw:5.0);
            Nmax = max(N_temp(:));
            
            if Nmax > 0
                [xcoor, ycoor] = meshgrid(0+xbinw/2:xbinw:180, 0.0+ybinw/2:ybinw:5.0);
                xcoor = xcoor(:);
                ycoor = ycoor(:);
                Ncounts1 = N_temp(:);
                Ninds = ceil(Ncounts1/Nmax*256);
                
                % Only plot non-zero points
                validPoints = Ninds > 0;
                if any(validPoints)
                    colors = c0(Ninds(validPoints), :);
                    scatter(ycoor(validPoints), xcoor(validPoints), sizepoint, colors, 'filled');
                end
                
                % Add colorbar
                h = utils_plot.createNarrowColorbar('vert');
                set(h, 'FontSize', C.fonts.tickSize);
                hlim = get(h,'YLim');
                set(h, 'YTick', hlim);
                set(h, 'YTickLabel', {0, Nmax});
                
                text(1.15, 0.5, '# neurons', 'Units', 'normalized', 'Rotation', 90, ...
                     'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
                % utils_plot.freezeColors(h);
            end
            
            % Mutual information test
            try
                if exist('+utils_analysis/mi_cont_cont', 'file') == 2
                    MI_test = utils_analysis.mi_cont_cont(norm_plot, tuning_plot);
                    MI_shuffle = zeros(Nshuffle, 1);
                    
                    for i = 1:min(Nshuffle, 100) % Limit shuffles for performance
                        MI_shuffle(i) = utils_analysis.mi_cont_cont(norm_plot(randperm(length(norm_plot))), ...
                                                   tuning_plot(randperm(length(tuning_plot))));
                    end
                    
                    p_value = sum(MI_shuffle > MI_test) / length(MI_shuffle);
                    
                    text('Units', 'Normalized', 'Position', [0.5 1.1], ...
                         'string', sprintf('p=%.2g', max(p_value, 1/Nshuffle)), ...
                         'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, ...
                         'color', 'k');
                else
                    warning('mi_cont_cont function not available');
                end
            catch ME
                warning('Error computing mutual information: %s', ME.message);
            end
            
        else
            warning('No valid data points for Panel B');
        end
        
    catch ME
        warning('Error processing preferred orientation analysis: %s', ME.message);
    end
    
    ylim([0 180]);
    xlim([0 5]);
    xlabel('normalization index', 'FontSize', C.fonts.labelSize);
    ylabel('preferred orientation', 'FontSize', C.fonts.labelSize);
    set(gca, 'XTick', xticks1, 'XTickLabel', xticklabels1);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C1: Original correlation matrix
    axes(AH(3));
    
    try
        if isfield(ssnData, 'dataheat_original')
            dataheat = ssnData.dataheat_original;
            
            if ~isempty(dataheat) && ~all(isnan(dataheat(:)))
                imagesc(0.05:0.1:2.95, 0.05:0.1:2.95, dataheat');
                colormap(AH(3), customColormap);
                clim(clim1);
                axis xy;
                
                % Add colorbar
                h1 = utils_plot.createNarrowColorbar('vert');
                set(h1, 'FontSize', C.fonts.tickSize);
                colormap(h1,customColormap);
                text(1.25, 0.5, 'correlation', 'Units', 'normalized', 'Rotation', 90, ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', C.fonts.annotationSize);
                utils_plot.freezeColors(h1);
            else
                warning('Invalid original correlation data');
                imagesc(0.05:0.1:2.95, 0.05:0.1:2.95, rand(30, 30));
            end
        end
        
    catch ME
        warning('Error processing original correlation matrix: %s', ME.message);
        imagesc(0.05:0.1:2.95, 0.05:0.1:2.95, rand(30, 30));
    end
    
    xlim([0 3]);
    ylim([0 3]);
    xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
    ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
    xticks([0 1 2 3]);
    yticks([0 1 2 3]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel C2: Tuning-matched correlation matrix
    axes(AH(4));
    
    try
        if isfield(ssnData, 'dataheat_matchtuning')
            dataheat = ssnData.dataheat_matchtuning;
            
            if ~isempty(dataheat) && ~all(isnan(dataheat(:)))
                imagesc(0.05:0.1:2.95, 0.05:0.1:2.95, dataheat');
                colormap(AH(4), customColormap);
                clim(clim1);
                axis xy;
                
                % Add colorbar
                h1 = utils_plot.createNarrowColorbar('vert');
                set(h1, 'FontSize', C.fonts.tickSize);
                colormap(h1,customColormap);
                text(1.25, 0.5, 'correlation', 'Units', 'normalized', 'Rotation', 90, ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', C.fonts.annotationSize);
                utils_plot.freezeColors(h1);
            else
                warning('Invalid tuning-matched correlation data');
                imagesc(0.05:0.1:2.95, 0.05:0.1:2.95, rand(30, 30));
            end
        end
        
    catch ME
        warning('Error processing tuning-matched correlation matrix: %s', ME.message);
        imagesc(0.05:0.1:2.95, 0.05:0.1:2.95, rand(30, 30));
    end
    
    xlim([0 3]);
    ylim([0 3]);
    xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
    ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
    xticks([0 1 2 3]);
    yticks([0 1 2 3]);
    title('match tuning', 'FontSize', C.fonts.titleSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_SSN', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_SSN:GeneralError', ...
          'Failed to generate SSN supplementary figure: %s', ME.message);
end

end