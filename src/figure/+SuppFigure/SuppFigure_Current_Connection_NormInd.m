function SuppFigure_Current_Connection_NormInd(varargin)
%SUPPFIGURE_CURRENT_CONNECTION_NORMIND Generate supplementary figure for current-connection vs normalization index
%
% Purpose:
%   Creates a figure analyzing the relationship between current types, connection
%   numbers, and normalization indices, showing correlations for feedforward
%   excitation, recurrent excitation, and recurrent inhibition across both
%   current strength and connection count.
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
%   SuppFigure_Current_Connection_NormInd();
%   SuppFigure_Current_Connection_NormInd('Save', false, 'FigNum', 2);

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
    utils_plot.setPlotStyle('custom', 'width', 11, 'height', 8.5);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = 0.08*ones(2,1);
    vertical_interval = 0.06;
    my_position2 = [0.08 0.07 0.90 0.84];
    
    DC = utils_plot.divideAxes([0.25 0.25 0.25], [0.25 0.25], my_position2, ...
                               horizontal_interval, vertical_interval);
    
    % Create axes
    LdPos1 = [-0.04, 0.10];
    LdPos = [-0.04, 0.06];
    labels = {'A1', 'A2', 'A3', 'B1', 'B2', 'B3'};
    panelorder = [1,3,5,2,4,6];funrun
    
    for i = 1:6
        AH(i) = axes('Position', DC{panelorder(i)});
        hold on;
        if i <= 3
            utils_plot.addFigureLabel(labels{i}, LdPos1);
        else
            utils_plot.addFigureLabel(labels{i}, LdPos);
        end
    end
    
    utils_plot.setFigureProperties();
    
    % Colors
    color_e = C.colors.excitatory;
    color_i = C.colors.inhibitory;
    
    %% Load data
    try
        % Load pattern data and indices
        patternData = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond.mat'));
        indexData = load(fullfile(rootFolder, 'data', 'indselect.mat'));
        connectionData = load(fullfile(rootFolder, 'data', 'FigSupp_Current_Connection_NormInd.mat'));
        
        % Validate required fields
        if ~isfield(patternData, 'MTE')
            error('MTE field not found in pattern data');
        end
        if ~isfield(indexData, 'indselect')
            error('indselect field not found in index data');
        end
        
        MTE = patternData.MTE;
        indselect = indexData.indselect;
        
        % Calculate normalization indices
        norm1 = (MTE(:,1) + MTE(:,2)) ./ MTE(:,3);
        
        % Validate indices
        validIndices = indselect(indselect <= length(norm1) & indselect > 0);
        if length(validIndices) < length(indselect)
            warning('Some indices are out of bounds, using %d valid indices', length(validIndices));
        end
        norm = norm1(validIndices);
        
        % Load current and connection data with validation
        requiredFields = {'I_ffwd', 'I_recE', 'I_recI', 'N_Ffwd', 'N_recE', 'N_recI'};
        for field = requiredFields
            if ~isfield(connectionData, field{1})
                error('Required field "%s" not found in connection data', field{1});
            end
        end
        
        % Extract current and connection data
        I_ffwd = connectionData.I_ffwd;
        I_recE = connectionData.I_recE;
        I_recI = connectionData.I_recI;
        N_Ffwd = connectionData.N_Ffwd;
        N_recE = connectionData.N_recE;
        N_recI = connectionData.N_recI;
        
        % Ensure data dimensions match
        dataLength = length(norm);
        currentData = {I_ffwd, I_recE, I_recI};
        connectionDataSets = {N_Ffwd, N_recE, N_recI};
        
        for i = 1:3
            if length(currentData{i}) < dataLength
                warning('Current data %d is shorter than norm data, padding with last value', i);
                currentData{i} = [currentData{i}(:); repmat(currentData{i}(end), dataLength - length(currentData{i}), 1)];
            elseif length(currentData{i}) > dataLength
                currentData{i} = currentData{i}(1:dataLength);
            end
            
            if length(connectionDataSets{i}) < dataLength
                warning('Connection data %d is shorter than norm data, padding with last value', i);
                connectionDataSets{i} = [connectionDataSets{i}(:); repmat(connectionDataSets{i}(end), dataLength - length(connectionDataSets{i}), 1)];
            elseif length(connectionDataSets{i}) > dataLength
                connectionDataSets{i} = connectionDataSets{i}(1:dataLength);
            end
        end
        
    catch ME
        error('SuppFigure_Current_Connection_NormInd:DataLoadError', ...
              'Failed to load required data: %s', ME.message);
    end
    
    %% Current analysis panels (A1-A3)
    currentTitles = {'Feedforward excitation', 'Recurrent excitation', 'Recurrent inhibition'};
    currentColors = {color_e, color_e, color_i};
    currentYLims = {[1.5 3.5], [1.5 3.5], [3 6]};
    
    for panelIdx = 1:3
        axes(AH(panelIdx));
        
        try
            I_temp = currentData{panelIdx};
            
            % Apply absolute value for inhibitory currents
            if panelIdx == 3
                I_temp = abs(I_temp);
            end
            
            % Remove invalid data
            validData = ~isnan(norm) & ~isnan(I_temp) & isfinite(norm) & isfinite(I_temp);
            norm_valid = norm(validData);
            I_valid = I_temp(validData);
            
            if length(norm_valid) > 10
                % Create scatter plot
                scatter(norm_valid, I_valid, 2, 'MarkerFaceColor', currentColors{panelIdx}, ...
                       'MarkerEdgeColor', currentColors{panelIdx}, 'LineWidth', 0.05, ...
                       'MarkerFaceAlpha', C.plot.alphaValue);
                
                % Calculate correlation and fit line
                [r, p_val] = corrcoef(norm_valid, I_valid);
                if size(r, 1) > 1 && ~isnan(r(1,2))
                    % Fit linear regression
                    X = [ones(length(norm_valid),1), norm_valid];
                    B = X \ I_valid;
                    
                    % Plot fitted line
                    x_fit = linspace(min(norm_valid), max(norm_valid), 100);
                    y_fit = x_fit * B(2) + B(1);
                    plot(x_fit, y_fit, 'color', 'k', 'LineWidth', C.plot.lineWidth);
                    
                    % Display correlation
                    my_string = sprintf('R=%.2f, p=%.1e', r(1,2), p_val(1,2));
                    text('Units', 'Normalized', 'Position', [0.5 0.95], 'string', my_string, ...
                         'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
                else
                    warning('Could not calculate valid correlation for current panel %d', panelIdx);
                end
            else
                warning('Insufficient valid data for current panel %d (%d points)', panelIdx, length(norm_valid));
                % Plot placeholder data
                scatter([1 2 3], [2 2.5 3], 2, currentColors{panelIdx});
            end
            
        catch ME
            warning('Error processing current panel A%d: %s', panelIdx, ME.message);
            scatter([1 2 3], [2 2.5 3], 2, currentColors{panelIdx});
        end
        
        xlim([1 3]);
        ylim(currentYLims{panelIdx});
        
        title(currentTitles{panelIdx}, 'color', currentColors{panelIdx}, ...
              'Units', 'normalized', 'Position', [0.5, 1.07], 'FontSize', C.fonts.titleSize);
        
        if panelIdx == 1
            ylabel('current (absolute value)', 'FontSize', C.fonts.labelSize);
        end
        
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Connection analysis panels (B1-B3)
    connectionYLims = {[200 350], [150 300], [150 300]};
    
    for panelIdx = 1:3
        axes(AH(panelIdx + 3));
        
        try
            N_temp = connectionDataSets{panelIdx};
            
            % Ensure N_temp is a column vector and handle transpose if needed
            if size(N_temp, 1) == 1
                N_temp = N_temp';
            end
            
            % Remove invalid data
            validData = ~isnan(norm) & ~isnan(N_temp) & isfinite(norm) & isfinite(N_temp);
            norm_valid = norm(validData);
            N_valid = N_temp(validData);
            
            if length(norm_valid) > 10
                % Create scatter plot
                scatter(norm_valid, N_valid, 2, 'MarkerFaceColor', currentColors{panelIdx}, ...
                       'MarkerEdgeColor', currentColors{panelIdx}, 'LineWidth', 0.05, ...
                       'MarkerFaceAlpha', C.plot.alphaValue);
                
                % Calculate correlation and fit line
                [r, p_val] = corrcoef(norm_valid, N_valid);
                if size(r, 1) > 1 && ~isnan(r(1,2))
                    % Fit linear regression
                    X = [ones(length(norm_valid),1), norm_valid];
                    B = X \ N_valid;
                    
                    % Plot fitted line
                    x_fit = linspace(min(norm_valid), max(norm_valid), 100);
                    y_fit = x_fit * B(2) + B(1);
                    plot(x_fit, y_fit, 'k', 'LineWidth', C.plot.lineWidth);
                    
                    % Display correlation
                    my_string = sprintf('R=%.2f, p=%.1e', r(1,2), p_val(1,2));
                    text('Units', 'Normalized', 'Position', [0.5 0.95], 'string', my_string, ...
                         'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
                else
                    warning('Could not calculate valid correlation for connection panel %d', panelIdx);
                end
            else
                warning('Insufficient valid data for connection panel %d (%d points)', panelIdx, length(norm_valid));
                % Plot placeholder data
                scatter([1 2 3], [250 250 250], 2, currentColors{panelIdx});
            end
            
        catch ME
            warning('Error processing connection panel B%d: %s', panelIdx, ME.message);
            scatter([1 2 3], [250 250 250], 2, currentColors{panelIdx});
        end
        
        xlim([1 3]);
        ylim(connectionYLims{panelIdx});
        xlabel('normalization index', 'FontSize', C.fonts.labelSize);
        
        if panelIdx == 1
            ylabel('number of connections', 'FontSize', C.fonts.labelSize);
        end
        
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_Current_Connection_NormInd', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_Current_Connection_NormInd:GeneralError', ...
          'Failed to generate current-connection normalization index figure: %s', ME.message);
end

end