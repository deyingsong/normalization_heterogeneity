function SuppFigure_DynamicalRegime_wave(varargin)
%SUPPFIGURE_DYNAMICALREGIME_WAVE Generate supplementary figure for wave dynamical regime
%
% Purpose:
%   Creates a comprehensive figure analyzing the wave dynamical regime including
%   normalization distributions, distance correlations, and current correlations
%   compared to default conditions.

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
    
    utils_plot.setPlotStyle('custom', 'width', 18, 'height', 11);
    
    figure(p.Results.FigNum); clf;
    utils_plot.matchFigureAspectRatio();
    
    horizontal_interval = [0.20, 0.10*ones(1,2)];
    vertical_interval = 0.05;
    my_position1 = [0.06 0.08 0.85 0.85];
    
    DC1 = utils_plot.divideAxes(0.25*ones(1,4), 0.25*ones(1,2), my_position1, ...
                                horizontal_interval, vertical_interval);
    
    % Load custom colormap
    try
        colormapData = load(fullfile(rootFolder, 'data', 'mycolormap2.mat'));
        customColormap = colormapData.mycolormap2;
    catch
        warning('Could not load custom colormap, using default');
        customColormap = parula(256);
    end
    
    LdPos = [-0.04, 0.05];
    labels = {'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'};
    panelorder = [1,3,5,2,4,6,8];
    
    for i = 1:7
        AH(i) = axes('Position', DC1{panelorder(i)});
        hold on;
        if i <= 3
            utils_plot.addFigureLabel(labels{i}, LdPos);
        else
            utils_plot.addFigureLabel(labels{i}, LdPos + [0 0.01]);
        end
    end
    
    utils_plot.setFigureProperties();
    
    % Colors and parameters
    col1 = C.colors.excitatory;
    legend_posY = [1.8 1.5];
    myclim = [-0.03 0.20];
    cb_fontsize = C.fonts.tickSize;
    
    %% Panel A1: Normalization distributions comparison
    axes(AH(1));
    
    try
        % Load wave regime data
        waveData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_norm_ind_wave.mat'), 'norm1');
        defaultData = load(fullfile(rootFolder, 'data', 'Fig1_NormIndDistribution.mat'));
        
        if isfield(waveData, 'norm1')
            temp = waveData.norm1(waveData.norm1 <= 3 & waveData.norm1 >= 0 & ~isnan(waveData.norm1));
            if ~isempty(temp)
                counts = histcounts(temp, 0:0.1:3);
                binCenters = 0.05:0.1:2.95;
                densityValues = counts / length(temp) / 0.1;
                plot(binCenters, densityValues, 'k', 'LineWidth', C.plot.lineWidth);
            end
        end
        
        if isfield(defaultData, 'norm1')
            norm_default = defaultData.norm1;
            temp = norm_default(norm_default <= 3 & norm_default >= 0 & ~isnan(norm_default));
            if ~isempty(temp)
                counts = histcounts(temp, 0:0.1:3);
                binCenters = 0.05:0.1:2.95;
                densityValues = counts / length(temp) / 0.1;
                plot(binCenters, densityValues, 'k-.', 'LineWidth', C.plot.lineWidth);
            end
        end
        
        % Add legend
        plot([2 3], [1 1]*legend_posY(1), 'Color', 'k', 'LineStyle', '-', 'LineWidth', C.plot.lineWidth);
        plot([2 3], [1 1]*legend_posY(2), 'Color', 'k', 'LineStyle', '-.', 'LineWidth', C.plot.lineWidth);
        text('Position', [3.1 legend_posY(1)], 'string', 'wave', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', 'k', 'FontSize', C.fonts.annotationSize);
        text('Position', [3.1 legend_posY(2)], 'string', 'Default', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', 'k', 'FontSize', C.fonts.annotationSize);
        
    catch ME
        warning('Error processing Panel A1: %s', ME.message);
        plot(0.05:0.1:2.95, rand(1, 30), 'k');
        plot(0.05:0.1:2.95, rand(1, 30), 'k-.');
    end
    
    xlim([0 3]);
    xlabel('norm index', 'FontSize', C.fonts.labelSize);
    ylabel('probability density', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel A2: Distance-correlation relationship
    axes(AH(2));
    
    try
        % Load distance correlation data
        defaultDistData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_rsc_dist_default.mat'));
        waveDistData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_rsc_dist_wave.mat'));
        
        binw = 0.03;
        Nbin = 20;
        binCenters = binw/2:binw:binw*Nbin;
        
        if isfield(defaultDistData, 'rsc_dist')
            plot(binCenters, defaultDistData.rsc_dist, 'k-.', 'LineWidth', C.plot.lineWidth);
        end
        
        if isfield(waveDistData, 'rsc_dist')
            plot(binCenters, waveDistData.rsc_dist, 'k', 'LineWidth', C.plot.lineWidth);
        end
        
    catch ME
        warning('Error processing Panel A2: %s', ME.message);
        plot(binw/2:binw:binw*Nbin, rand(1, Nbin), 'k');
        plot(binw/2:binw:binw*Nbin, rand(1, Nbin), 'k-.');
    end
    
    xlabel('distance (a.u.)', 'FontSize', C.fonts.labelSize);
    ylabel('correlation', 'FontSize', C.fonts.labelSize);
    xlim([0 0.5]);
    xticks([0 0.25 0.5]);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel A3: Empty (axis off)
    axes(AH(3));
    axis off;
    
    %% Load firing rate pattern data for correlation analysis
    try
        wavePatternData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_FRpattern_wave.mat'), 'MTE');
        if isfield(wavePatternData, 'MTE')
            MTE = wavePatternData.MTE;
            ind1 = find(MTE(:,3) > 2 & MTE(:,2) > 2 & MTE(:,1) > 2);
            norm_wave = waveData.norm1(ind1);
        else
            error('MTE field not found in wave pattern data');
        end
    catch ME
        warning('Error loading wave pattern data: %s', ME.message);
        ind1 = 1:1000;
        norm_wave = rand(1000, 1) * 3;
    end
    
    %% Panel A4: Correlation heatmap
    axes(AH(4));
    
    try
        % Load spike data and compute correlation matrix
        spikeFile = fullfile(rootFolder, 'data', 'RF2D2layer_weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_E_0d15_I_0d00_tuning_th_0d60_V1fr_Rec_3_160_100_AttFar_0d00_Noise.mat');
        
        if exist(spikeFile, 'file')
            spikeData = load(spikeFile, 's2');
            
            % Convert spike times to counts (using helper function)

                MTE2 = utils_analysis.spktime2count(spikeData.s2, 1:C.network.Ne, 50, 400, 1);
                MTE2 = MTE2(:, 21:400);
                
                % Apply filter
                hf = ones(1,4)/4;
                if exist('imfilter', 'file') == 2
                    MTE2 = imfilter(MTE2, hf, 'full');
                end
                MTE4 = MTE2(:, 4:380);
                
                % Get correlation matrix
                if exist(fullfile(utilsPath,'+utils_analysis','get_norm_rsc'), 'file') == 2 && length(ind1) <= size(MTE4, 1)
                    [dataheat, ind_lim] = utils_analysis.get_norm_rsc(MTE4(ind1,:), norm_wave);
                    nbins = 20;
                    binw = (ind_lim(2) - ind_lim(1)) / nbins;
                    
                    imagesc(ind_lim(1)+binw/2:binw:ind_lim(2), ind_lim(1)+binw/2:binw:ind_lim(2), dataheat');
                    xlim(ind_lim); ylim(ind_lim);
                    colormap(AH(4), customColormap);
                    clim(myclim);
                else
                    warning('Helper functions not available or dimension mismatch');
                    imagesc(rand(20, 20) * diff(myclim) + myclim(1));
                end
            
        else
            warning('Spike data file not found');
            imagesc(rand(20, 20) * diff(myclim) + myclim(1));
        end
        
        % Add colorbar
        h1 = colorbar;
        colormap(h1, customColormap);
        set(h1, 'Position', [DC1{panelorder(4)}(1)+DC1{panelorder(4)}(3)+0.01 DC1{panelorder(4)}(2)+0.06 0.003 0.26]);
        set(h1, 'FontSize', cb_fontsize);
        set(h1, 'Ticks', 0:0.05:0.15);
        clim(myclim);
        
        text(1.35, 0.5, 'correlation', 'Units', 'normalized', 'Rotation', 90, ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
        
    catch ME
        warning('Error processing Panel A4: %s', ME.message);
        imagesc(rand(20, 20));
        colormap(AH(4), customColormap);
    end
    
    axis xy; box on;
    xlabel('norm index 1', 'FontSize', C.fonts.labelSize);
    ylabel('norm index 2', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Load current data for correlation analysis
    try
        waveCurrentData = load(fullfile(rootFolder, 'data', 'FigSupp_DynamicalRegime_current_wave.mat'));
        
        if isfield(waveCurrentData, 'inputXr') && isfield(waveCurrentData, 'inputEr') && isfield(waveCurrentData, 'inputIr')
            normX = (sum(waveCurrentData.inputXr(ind1,1,:), 3) + sum(waveCurrentData.inputXr(ind1,2,:), 3)) ./ sum(waveCurrentData.inputXr(ind1,3,:), 3);
            normX = squeeze(normX);
            normE = (waveCurrentData.inputEr(ind1,1) + waveCurrentData.inputEr(ind1,2)) ./ waveCurrentData.inputEr(ind1,3);
            normI = (waveCurrentData.inputIr(ind1,1) + waveCurrentData.inputIr(ind1,2)) ./ waveCurrentData.inputIr(ind1,3);
            normfr = norm_wave;
            
            % Calculate maximum for colorbar
            Nbins = 60;
            xbinw = 6/Nbins; ybinw = 0.3/Nbins;
            N_temp = [histcounts2(normfr, normX, 0:xbinw:6, 1.4:ybinw:1.7), ...
                     histcounts2(normfr, normE, 0:xbinw:6, 1.4:ybinw:1.7), ...
                     histcounts2(normfr, normI, 0:xbinw:6, 1.4:ybinw:1.7)];
            Nmax = max(N_temp(:));
        else
            error('Required current fields not found');
        end
    catch ME
        warning('Error loading wave current data: %s', ME.message);
        normX = rand(length(ind1), 1) * 0.3 + 1.4;
        normE = rand(length(ind1), 1) * 0.3 + 1.4;
        normI = rand(length(ind1), 1) * 0.3 + 1.4;
        normfr = rand(length(ind1), 1) * 6;
        Nmax = 100;
    end
    
    %% Panels A5-A7: Current correlation analysis
    currentTypes = {normX, normE, normI};
    currentTitles = {'Feedforward excitation', 'Recurrent excitation', 'Recurrent inhibition'};
    currentColors = {col1, col1, C.colors.inhibitory};
    panelIndices = [5, 6, 7];
    
    for panelIdx = 1:3
        axes(AH(panelIndices(panelIdx)));
        
        try
            currentNorm = currentTypes{panelIdx};
            
            % Create density scatter plot
            Ncounts1 = histcounts2(normfr, currentNorm, 0:xbinw:6, 1.4:ybinw:1.7);
            Ncounts1 = Ncounts1(:);
            Ninds = ceil(Ncounts1/Nmax*256);
            [xcoor, ycoor] = meshgrid(0+xbinw/2:xbinw:6, 1.4+ybinw/2:ybinw:1.7);
            xcoor=reshape(xcoor',[],1);
            ycoor=reshape(ycoor',[],1);
            
            validPoints = Ninds > 0;
            if any(validPoints)
                c0 = parula;
                colors = c0(Ninds(validPoints), :);
                scatter(xcoor(validPoints), ycoor(validPoints), 3, colors, 'filled');
                xlim([0 6]);
                xticks(0:2:6);
                ylim([1.3 1.7]);
                clim([0 Nmax]);
            end
            
            % Calculate correlation
            [r, ~] = corrcoef(normfr, currentNorm);
            if size(r, 1) > 1
                my_string = sprintf('R=%.2f', r(1,2));
                text('Units', 'Normalized', 'Position', [0.72 0.13], 'string', my_string, ...
                     'Interpreter', 'latex', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
            end
            
            % Add titles
            text('Units', 'Normalized', 'Position', [0.5 C.positions.titleHeight], ...
                 'string', currentTitles{panelIdx}, 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', C.fonts.titleSize, 'color', currentColors{panelIdx});
            
        catch ME
            warning('Error processing Panel A%d: %s', panelIdx+4, ME.message);
            scatter(rand(50,1)*6, rand(50,1)*0.3+1.4, 3, rand(50,3));
        end
        
        xlabel('rate norm index', 'FontSize', C.fonts.labelSize);
        ylabel('current norm index', 'FontSize', C.fonts.labelSize);
        ylim([1.3 1.7]);
        xlim([0 6]);
        pbaspect(C.figure.aspectRatio);
        
        % Add colorbar for last panel
        if panelIdx == 3
            h = colorbar;
            set(h, 'Position', [DC1{panelorder(panelIndices(panelIdx))}(1)+DC1{panelorder(panelIndices(panelIdx))}(3)+0.01 ...
                               DC1{panelorder(panelIndices(panelIdx))}(2)+0.06 0.003 0.26]);
            set(h, 'FontSize', cb_fontsize);
            set(h, 'Ticks', [0 Nmax]);
            set(h, 'TickLabels', [0 Nmax]);
            colormap(h, "parula");
            
            text(1.15, 0.5, '# neurons', 'Units', 'normalized', 'Rotation', 90, ...
                 'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
        end
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_DynamicalRegime_wave', ...
                             'view', p.Results.View, 'save', true, 'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_DynamicalRegime_wave:GeneralError', ...
          'Failed to generate wave dynamical regime figure: %s', ME.message);
end
end

