function SuppFigure_WeakCoupling_EIbalance(varargin)
%SUPPFIGURE_WEAKCOUPLING_EIBALANCE Generate supplementary figure for weak coupling E/I balance analysis
%
% Purpose:
%   Creates a comprehensive figure showing the relationship between coupling
%   strength and excitation-inhibition balance, including balance indices,
%   normalization distributions, and correlation heatmaps across different
%   coupling scales.
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
%   SuppFigure_WeakCoupling_EIbalance();
%   SuppFigure_WeakCoupling_EIbalance('Save', false, 'FigNum', 2);

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
    utils_plot.setPlotStyle('custom', 'width', C.figure.width + 1, 'height', C.figure.height + 4);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Define layout
    horizontal_interval = [0.09 0.09];
    vertical_interval = [0.12 0.12];
    my_position1 = [0.06 0.08 0.85 0.85];
    
    DC1 = utils_plot.divideAxes(0.25*ones(1,3), 0.25*ones(1,3), my_position1, ...
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
    
    % Create axes with labels
    LdPos = C.positions.legendOffset;
    labels = {'A', 'B', 'C1', 'C2', 'C3', 'C4', 'C5'};
    panelorder = [1,4,2,5,8,3,6];
    
    for i = 1:7
        AH(i) = axes('Position', DC1{panelorder(i)});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties();
    
    % Define parameters
    scaleJs = [0.1, 0.25, 0.5, 0.75, 1, 2];
    scaledown_ffwds = [0.34, 0.45, 0.63, 0.83, 1.0, 2.0];
    sigma_current = 6.8;
    att_cond = 'AttFar';
    inI = 0;
    sim_flag = 'CurrentNoise';
    
    % Color scheme
    color1 = cool(5);
    color2 = color1(1:4,:);
    color2(5,:) = [0 0 0];
    color2(6,:) = color1(5,:);
    mysz = 20;
    
    labels_scale = {'0.1', '0.25', '0.5', '0.75', '1(default)', '2'};
    myclims = [0.02, 0.04; 0.03, 0.06; 0.02, 0.05; 0, 0.04; 0, 0.01];
    
    %% Panel A: Balance index vs scale
    axes(AH(1));
    
    Jorder = [1,2,3,4,6];
    try
        mu1 = zeros(length(scaleJs), 1);
        
        for i = 1:length(scaleJs)
            scaleJ = 1/scaleJs(i);
            scaledown_ffwd = 1/scaledown_ffwds(i);
            
            % Adjust simulation flag based on scale
            current_sim_flag = sim_flag;
            current_sigma = sigma_current;
            if scaleJ <= 1 && strcmp(sim_flag, 'CurrentNoise')
                current_sim_flag = 'Noise';
                current_sigma = 0;
            end
            
            % Generate filename
            filename = sprintf('Mean_Current_weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_E_0d15_I_0d00_tuning_th_0d60_V1fr_Rec_3_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%.2f_%s', ...
                              1/scaleJ, 1/scaledown_ffwd, current_sigma, att_cond, inI, current_sim_flag);
            filename = strrep(filename, '.', 'd');
            
            try
                data1 = load(fullfile(rootFolder, 'data', [filename '.mat']), 'mu');
                if isfield(data1, 'mu') && ~isnan(data1.mu)
                    mu1(i) = data1.mu;
                else
                    warning('Invalid mu value for scale %d', i);
                    mu1(i) = NaN;
                end
            catch
                warning('Could not load data for scale index %d', i);
                mu1(i) = NaN;
            end
        end
        
        % Plot only valid data points
        validIndices = ~isnan(mu1);
        if any(validIndices)
            plot(scaleJs(validIndices), mu1(validIndices), 'k', 'LineWidth', C.plot.lineWidth);
            
            for i = 1:length(scaleJs)
                if validIndices(i)
                    scatter(scaleJs(i), mu1(i), mysz, color2(i,:));
                end
            end
        end
        
    catch ME
        warning('Error processing balance index data: %s', ME.message);
        plot(scaleJs, zeros(size(scaleJs)), 'k');
    end
    
    xlabel('scale', 'FontSize', C.fonts.labelSize);
    ylabel('balance index', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Normalization index distributions
    axes(AH(2));
    
    try
        for i = 1:length(scaleJs)
            scaleJ = 1/scaleJs(i);
            scaledown_ffwd = 1/scaledown_ffwds(i);
            
            % Adjust simulation parameters
            current_sim_flag = sim_flag;
            current_sigma = sigma_current;
            if scaleJ <= 1 && strcmp(sim_flag, 'CurrentNoise')
                current_sim_flag = 'Noise';
                current_sigma = 0;
            end
            
            filename = sprintf('Norm_ind_weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_E_0d15_I_0d00_tuning_th_0d60_V1fr_Rec_3_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%.2f_%s', ...
                              1/scaleJ, 1/scaledown_ffwd, current_sigma, att_cond, inI, current_sim_flag);
            filename = strrep(filename, '.', 'd');
            
            try
                data1 = load(fullfile(rootFolder, 'data', [filename '.mat']), 'pdf1');
                if isfield(data1, 'pdf1') && ~isempty(data1.pdf1)
                    xvals = 0.05:0.1:3;
                    if length(data1.pdf1) == length(xvals)
                        plot(xvals, data1.pdf1, 'Color', color2(i,:), 'LineWidth', C.plot.lineWidth);
                    end
                end
            catch
                warning('Could not load normalization distribution for scale %d', i);
            end
            
            % Add legend text
            text('Units', 'Normalized', 'Position', [1.3 0.95-i*0.12], ...
                 'string', labels_scale{i}, 'Interpreter', 'latex', ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                 'FontSize', C.fonts.annotationSize, 'color', color2(i,:));
        end
        
        % Legend title
        text('Units', 'Normalized', 'Position', [1.3 0.95], ...
             'string', 'Scale of recurrent connection', 'Interpreter', 'latex', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', C.fonts.annotationSize, 'color', 'k');
        
    catch ME
        warning('Error processing normalization distributions: %s', ME.message);
    end
    
    xlabel('normalization index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panels C1-C5: Correlation heatmaps
    iAs = [3, 4, 5, 6, 7];
    inds = [1, 2, 3, 4, 6]; % Skip index 5
    
    for panelIdx = 1:5
        axes(AH(iAs(panelIdx)));
        
        i = inds(panelIdx);
        scaleJ = 1/scaleJs(i);
        scaledown_ffwd = 1/scaledown_ffwds(i);
        
        % Adjust simulation parameters
        current_sim_flag = sim_flag;
        current_sigma = sigma_current;
        if scaleJ <= 1 && strcmp(sim_flag, 'CurrentNoise')
            current_sim_flag = 'Noise';
            current_sigma = 0;
        end
        
        filename = sprintf('Norm_rsc_preset_indlim_weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_E_0d15_I_0d00_tuning_th_0d60_V1fr_Rec_3_scaleJ_%.2f_scaleFFwd_%.2f_sigmaCurrent_%.1f_%s_%.2f_%s', ...
                          1/scaleJ, 1/scaledown_ffwd, current_sigma, att_cond, inI, current_sim_flag);
        filename = strrep(filename, '.', 'd');
        
        try
            data1 = load(fullfile(rootFolder, 'data', [filename '.mat']));
            
            if isfield(data1, 'dataheat') && isfield(data1, 'ind_lim')
                lower_lim = data1.ind_lim(1);
                upper_lim = data1.ind_lim(2);
                nbins = 20;
                binw = (upper_lim - lower_lim) / nbins;
                
                imagesc(lower_lim + binw/2:binw:upper_lim, ...
                       lower_lim + binw/2:binw:upper_lim, data1.dataheat');
                axis xy;
                
                % Set color limits safely
                clim_idx = min(panelIdx, size(myclims, 1));
                clim(myclims(clim_idx, :));
                colormap(AH(iAs(panelIdx)), customColormap);
                
                xlim([lower_lim upper_lim]);
                ylim([lower_lim upper_lim]);
                
                % Set appropriate ticks
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
                colormap(h1, customColormap);
                utils_plot.freezeColors(h1);
                
                % Add scale label
                text('Units', 'Normalized', 'Position', [0.5 1.12], ...
                     'string', labels_scale{i}, 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', C.fonts.titleSize, ...
                     'color', color2(i,:));
                
            else
                warning('Invalid data structure for panel C%d', panelIdx);
                imagesc(rand(20, 20)); % Placeholder
            end
            
        catch
            warning('Could not load correlation data for panel C%d', panelIdx);
            imagesc(rand(20, 20)); % Placeholder
        end
        
        xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
        if panelIdx == 1 || panelIdx == 4
            ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
        end
        
        if panelIdx == 3 || panelIdx == 5
            text('Units', 'Normalized', 'Position', [1.35 0.5], ...
                 'string', 'correlation', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, ...
                 'Rotation', 90, 'color', 'k');
        end
        
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_WeakCoupling_EIbalance', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_WeakCoupling_EIbalance:GeneralError', ...
          'Failed to generate supplementary figure: %s', ME.message);
end

end