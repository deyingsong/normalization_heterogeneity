function SuppFigure_VarOriPair(varargin)
%SUPPFIGURE_VARORIPAR Generate supplementary figure for variable orientation pair analysis
%
% Purpose:
%   Creates a figure showing the relationship between normalization indices
%   for different orientation pairs, comparing default model with superimposed
%   conditions and analyzing correlation coefficients across orientation differences.
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
%   SuppFigure_VarOriPair();
%   SuppFigure_VarOriPair('Save', false, 'FigNum', 2);

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
    utils_plot.setPlotStyle('custom', 'width', 12, 'height', 5);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = [0.18 0.08];
    vertical_interval = 0.08;
    my_position1 = [0.07 0.08 0.85 0.85];
    
    DC1 = utils_plot.divideAxes(0.25*ones(1,3), 0.5, my_position1, horizontal_interval, []);
    
    % Create axes
    LdPos = C.positions.legendOffset + [0 0.11];
    labels = {'A1', 'A2', 'B'};
    
    for i = 1:3
        AH(i) = axes('Position', DC1{i});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Parameters
    oris = [5,95; 10,100; 15,105; 20,110; 25,115; 30,120; 35,125; 40,130; 45,135];
    weight_name = 'weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_E_0d15_I_0d00_tuning_th_0d60';
    
    %% Panel A1: Default model correlation
    axes(AH(1));
    
    try
        % Load pattern data and normalization data
        patternData = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond.mat'));
        normData = load(fullfile(rootFolder, 'data', 'Norm_Current_I0_weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_E_0d15_I_0d00_tuning_th_0d60_V1fr_45_135_Rec.mat'), 'norm1');
        
        if ~isfield(patternData, 'MTE') || ~isfield(normData, 'norm1')
            error('Required fields not found in data files');
        end
        
        MTE = patternData.MTE;
        norm1_45_135 = normData.norm1;
        
        % Filter neurons based on firing rate criteria
        ind1 = find(abs(MTE(:,3) - mean(MTE(:,3))) <= std(MTE(:,3)) & ...
                   abs(MTE(:,1) - mean(MTE(:,1))) <= std(MTE(:,1)) & ...
                   abs(MTE(:,2) - mean(MTE(:,2))) <= std(MTE(:,2)));
        
        norm1_default = (MTE(:,1) + MTE(:,2)) ./ MTE(:,3);
        
        if length(ind1) > length(norm1_45_135)
            ind1 = ind1(1:length(norm1_45_135));
        end
        
        % Create density scatter plot
        Nbins = 60;
        xbinw = 3/Nbins;
        ybinw = 3/Nbins;
        sizepoint = 3;
        c0 = parula;
        
        N_temp = histcounts2(norm1_default(ind1), norm1_45_135(ind1), 0:xbinw:3, 0.0:ybinw:3.0);
        Nmax = max(N_temp(:));
        
        if Nmax > 0
            Ncounts1 = N_temp(:);
            Ninds = ceil(Ncounts1/Nmax*256);
            [xcoor,ycoor] = meshgrid(0+xbinw/2:xbinw:3.0, 0.0+ybinw/2:ybinw:3.0);
            xcoor = xcoor(:);
            ycoor = ycoor(:);
            
            validPoints = Ninds > 0;
            if any(validPoints)
                colors = c0(Ninds(validPoints), :);
                scatter(xcoor(validPoints), ycoor(validPoints), sizepoint, colors, 'filled');
            end
        end
        
        % Add diagonal line and formatting
        plot([0 3], [0 3], 'k:', 'LineWidth', C.plot.lineWidth);
        xlim([0 3]);
        ylim([0 3]);
        
        % Add colorbar
        adjust1 = [0.02 0.19 0 -0.38];
        h = utils_plot.createNarrowColorbar('vert');
        set(h, 'FontSize', C.fonts.tickSize);
        set(h, 'YTick', [0 1]);
        set(h, 'YTickLabel', [0 Nmax]);
        
        text(1.25, 0.5, '# neurons', 'Units', 'normalized', 'Rotation', 90, ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
        
        % Calculate and display correlation
        [r, ~] = corrcoef(norm1_default(ind1), norm1_45_135(ind1));
        text(0.5, 1.3, 'Default model', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.titleSize, 'color', 'k');
        text(0.5, 1.15, sprintf('R=%.2f', r(1,2)), 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.titleSize, 'color', 'k');
        
    catch ME
        warning('Error processing Panel A1: %s', ME.message);
        plot([0 3], [0 3], 'k');
        text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center');
    end
    
    xlabel('ori1=180°, ori2=90°', 'FontSize', C.fonts.labelSize);
    ylabel('ori1=45°, ori2=135°', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel A2: Correlation across orientation pairs
    axes(AH(2));
    
    try
        r1 = zeros(size(oris, 1), 1);
        ticklabels = cell(size(oris, 1), 1);
        
        for i = 1:size(oris, 1)
            V1_name = ['V1fr_' num2str(oris(i,1)) '_' num2str(oris(i,2)) '_Rec'];
            stats_name = sprintf('Norm_Current_I0_%s_%s', weight_name, V1_name);
            
            try
                oriData = load(fullfile(rootFolder, 'data', [stats_name '.mat']), 'norm1');
                if isfield(oriData, 'norm1')
                    [r, ~] = corrcoef(norm1_default(ind1), oriData.norm1(ind1));
                    r1(i) = r(1,2);
                else
                    r1(i) = NaN;
                end
            catch
                warning('Could not load data for orientation pair %d', i);
                r1(i) = NaN;
            end
            
            ticklabels{i} = num2str(oris(i,1));
        end
        
        % Plot correlation coefficients
        validCorr = ~isnan(r1);
        if any(validCorr)
            plot(1:size(oris,1), r1, 'k', 'LineWidth', C.plot.lineWidth);
        end
        
        set(gca, 'XTick', 1:9, 'XTickLabel', ticklabels);
        xlim([1 9]);
        xlabel('ori1', 'FontSize', C.fonts.labelSize);
        ylabel('correlation coefficient', 'FontSize', C.fonts.labelSize);
        
    catch ME
        warning('Error processing Panel A2: %s', ME.message);
        plot(1:9, zeros(1,9), 'k');
    end
    
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Superimposed condition correlation
    axes(AH(3));
    
    try
        % Load superimposed data
        data1 = load(fullfile(rootFolder, 'data', 'FigSupp_Superimposed_45_135_NormIndDistribution.mat'));
        data2 = load(fullfile(rootFolder, 'data', 'FigSupp_Superimposed_NormIndDistribution.mat'));
        controlData = load(fullfile(rootFolder, 'data', 'FigSupp_SuperimposedControl.mat'), 'MTE');
        
        if ~isfield(data1, 'norm1') || ~isfield(data2, 'norm1') || ~isfield(controlData, 'MTE')
            error('Required fields not found in superimposed data');
        end
        
        % Filter neurons
        MTE_control = controlData.MTE;
        ind1_super = find(abs(MTE_control(:,3) - mean(MTE_control(:,3))) <= std(MTE_control(:,3)) & ...
                         abs(MTE_control(:,1) - mean(MTE_control(:,1))) <= std(MTE_control(:,1)) & ...
                         abs(MTE_control(:,2) - mean(MTE_control(:,2))) <= std(MTE_control(:,2)));
        
        % Create density scatter plot
        N_temp = histcounts2(data2.norm1(ind1_super), data1.norm1(ind1_super), 0:xbinw:3, 0.0:ybinw:3.0);
        Nmax = max(N_temp(:));
        
        if Nmax > 0
            Ncounts1 = N_temp(:);
            Ninds = ceil(Ncounts1/Nmax*256);
            [xcoor,ycoor] = meshgrid(0+xbinw/2:xbinw:3.0, 0.0+ybinw/2:ybinw:3.0);
            xcoor = xcoor(:);
            ycoor = ycoor(:);
            
            validPoints = Ninds > 0;
            if any(validPoints)
                colors = c0(Ninds(validPoints), :);
                scatter(xcoor(validPoints), ycoor(validPoints), sizepoint, colors, 'filled');
            end
        end
        
        % Add diagonal and formatting
        plot([0 3], [0 3], 'k:', 'LineWidth', C.plot.lineWidth);
        xlim([0 3]);
        ylim([0 3]);
        
        % Add colorbar
        h = utils_plot.createNarrowColorbar('vert');
        set(h, 'FontSize', C.fonts.tickSize);
        set(h, 'YTick', [0 1]);
        set(h, 'YTickLabel', [0 Nmax]);
        
        text(1.25, 0.5, '# neurons', 'Units', 'normalized', 'Rotation', 90, ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
        
        % Calculate and display correlation
        [r, ~] = corrcoef(data2.norm1(ind1_super), data1.norm1(ind1_super));
        text(0.5, 1.3, 'Superimposed', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.titleSize, 'color', 'k');
        text(0.5, 1.15, sprintf('R=%.2f', r(1,2)), 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', C.fonts.titleSize, 'color', 'k');
        
    catch ME
        warning('Error processing Panel B: %s', ME.message);
        plot([0 3], [0 3], 'k');
        text(0.5, 0.5, 'Data not available', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center');
    end
    
    xlabel('ori1=180°, ori2=90°', 'FontSize', C.fonts.labelSize);
    ylabel('ori1=45°, ori2=135°', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_VarOriPair', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_VarOriPair:GeneralError', ...
          'Failed to generate variable orientation pair supplementary figure: %s', ME.message);
end

end