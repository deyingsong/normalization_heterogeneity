function SuppFigure_FFmodel(varargin)
%SUPPFIGURE_FFMODEL Generate supplementary figure for feedforward model analysis
%
% Purpose:
%   Creates a comprehensive figure showing feedforward model analysis including
%   firing rate vs total current relationships, normalization index patterns,
%   and theoretical predictions across different parameter regimes.
%
% Inputs:
%   varargin - Name-value pairs:
%     'Save'        - logical, whether to save figure (default: true)
%     'View'        - logical, whether to display figure (default: true)
%     'FigNum'      - integer, figure number (default: 10)
%     'Recompute'   - logical, whether to recompute data (default: false)
%
% Outputs:
%   None (creates and optionally saves figure)
%
% Usage Example:
%   SuppFigure_FFmodel();
%   SuppFigure_FFmodel('Recompute', true, 'FigNum', 5);

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
    addParameter(p, 'FigNum', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Recompute', false, @islogical);
    parse(p, varargin{:});
    
    % Set plot style
    utils_plot.setPlotStyle('custom', 'width', 17, 'height', 17);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Define layout
    DC = utils_plot.divideAxes(4, 4, [0.05 0.1 0.92 0.85], .3, .4);
    DC2 = utils_plot.divideAxes([0.8 1 1], 4, [0.1 0.05 0.8 0.85], .4, .4);
    panelorder=[2,6,10,3,7,11,4,8,12];
    
    % Create main axes
    for i = 1:4
        AH(i) = axes('Position', DC{i*4-3});
        hold on;
    end
    
    % Create secondary axes
    for i = 1:numel(DC2)-3
        AH2(i+3) = axes('Position', DC2{panelorder(i)});
        hold on;
    end
    
    % Create inset
    Pos = DC{5};
    Pos2 = Pos;
    Pos2(3) = Pos(3) * 0.4;
    Pos2(4) = Pos(4) * 0.4;
    Pos2(2) = Pos(2) + Pos(4) * 0.6;
    Pos2(1) = Pos(1) + Pos(3) * 0.1;
    AH_inset = axes('Position', Pos2);
    hold on;
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    %% Load and validate data
    try
        % Load main data file
        dataFile = fullfile(rootFolder, 'data', 'FigSupp_rate_current.mat');
        if ~exist(dataFile, 'file')
            error('Required data file not found: %s', dataFile);
        end
        
        data = load(dataFile);
        
        % Validate required fields
        requiredFields = {'Itot', 'rate', 'Ix', 'norm1'};
        for field = requiredFields
            if ~isfield(data, field{1})
                error('Required field "%s" not found in data file', field{1});
            end
        end
        
        Itot = data.Itot;
        rate = data.rate;
        Ix = data.Ix;
        norm1 = data.norm1;
        
        % Validate data dimensions
        if size(Itot, 2) < 3 || size(rate, 2) < 3
            error('Data must have at least 3 columns');
        end
        
        if length(norm1) ~= size(rate, 1)
            error('Mismatch between norm1 and rate data lengths');
        end
        
    catch ME
        error('SuppFigure_FFmodel:DataLoadError', ...
              'Failed to load required data: %s', ME.message);
    end
    
    %% Panel A: I1+I2 vs total current relationship
    axes(AH(1));
    
    try
        % Plot data points
        scatter(Itot(:,3), Itot(:,1) + Itot(:,2), 3, 'c', 'filled', ...
               'MarkerFaceAlpha', C.plot.alphaValue);
        hold on;
        scatter(Ix(:,3), Ix(:,1) + Ix(:,2), 3, 'm', 'filled', ...
               'MarkerFaceAlpha', C.plot.alphaValue);
        
        % Add reference line
        xlims = xlim;
        plot(xlims, xlims * 1.5, 'k--', 'LineWidth', C.plot.lineWidth);
        
    catch ME
        warning('Error plotting Panel A: %s', ME.message);
        plot([0 1], [0 1], 'k'); % Placeholder
    end
    
    xlabel('I_{1+2}', 'FontSize', C.fonts.labelSize);
    ylabel('I_1 + I_2', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Panel B: Firing rate vs total current with fits
    axes(AH(2));
    
    try
        % Plot raw data
        scatter(Itot(:,3), rate(:,3), 1, 'k', 'filled', ...
               'MarkerFaceAlpha', C.plot.alphaValue);
        
        % Fit piecewise function
        rate1 = rate(rate(:,3) < 20, 3);
        Itot1 = Itot(rate(:,3) < 20, 3);
        rate2 = rate(rate(:,3) >= 20, 3);
        Itot2 = Itot(rate(:,3) >= 20, 3);
        
        if ~isempty(rate1) && ~isempty(Itot1) && length(rate1) > 3
            % Fit first piece (low rates)
            try
                myfittype1 = fittype('k.*(x-a).^n.*(x>a)', ...
                                   'dependent', {'y'}, 'independent', {'x'}, ...
                                   'coefficients', {'k', 'a', 'n'});
                
                % Use robust starting points
                startPoint = [0.01, min(Itot1), 2];
                yE_both1 = fit(Itot1, rate1, myfittype1, 'StartPoint', startPoint);
                
                x0 = (20/yE_both1.k)^(1/yE_both1.n) + yE_both1.a;
                
                % Plot first piece
                x1 = linspace(min(Itot1), x0, 101);
                plot(x1, feval(yE_both1, x1), 'r', 'LineWidth', C.plot.lineWidth);
                
            catch fitError
                warning('Could not fit first piece: %s', fitError.message);
            end
        end
        
        if ~isempty(rate2) && ~isempty(Itot2) && length(rate2) > 3 && exist('x0', 'var')
            % Fit second piece (high rates)
            try
                myfittype2 = fittype('k.*(x-(x0 - (20/k).^(1/n))).^n.*(x>0)', ...
                                   'dependent', {'y'}, 'independent', {'x'}, ...
                                   'problem', 'x0', 'coefficients', {'k', 'n'});
                
                yE_both2 = fit(Itot2, rate2, myfittype2, 'problem', x0, ...
                              'StartPoint', [1, 1]);
                
                % Plot second piece
                x2 = linspace(x0, max(Itot2), 101);
                plot(x2, feval(yE_both2, x2), 'r', 'LineWidth', C.plot.lineWidth);
                
            catch fitError
                warning('Could not fit second piece: %s', fitError.message);
            end
        end
        
    catch ME
        warning('Error processing Panel B: %s', ME.message);
        plot(Itot(:,3), rate(:,3), '.'); % Fallback plot
    end
    
    xlabel('I_{tot}', 'FontSize', C.fonts.labelSize);
    ylabel('rate (Hz)', 'FontSize', C.fonts.labelSize);
    pbaspect(C.figure.aspectRatio);
    
    %% Inset: Zoom of low current region
    axes(AH_inset);
    
    try
        scatter(Itot(:,3), rate(:,3), 1, 'k', 'filled');
        if exist('x1', 'var') && exist('yE_both1', 'var')
            plot(x1, feval(yE_both1, x1), 'r', 'LineWidth', C.plot.lineWidth);
        end
        axis([-0.5 0.5 0 10]);
        
    catch ME
        warning('Error creating inset: %s', ME.message);
    end
    
    %% Panels C and D: Normalization patterns by current offset
    try
        % Calculate current offset
        % p1 = polyfit(Itot(:,1) + Itot(:,2), Itot(:,3), 1);
        th = Itot(:,3) - (Itot(:,1) + Itot(:,2))/1.5;
        bins = quantile(th, 0:0.25:1);
        
        plotPanels = [4, 1]; % Q4, Q1
        kk = 3;
        
        for bbIdx = 1:length(plotPanels)
            bb = plotPanels(bbIdx);
            axes(AH(bbIdx+2));
            
            try
                ind1 = find((rate(:,3) > 1) & (th < bins(bb+1)) & (th >= bins(bb)));
                
                if ~isempty(ind1)
                    % Create scatter plot colored by normalization index
                    scatter(rate(ind1,1), rate(ind1,2), 5, norm1(ind1), 'filled');
                    axis([0 80 0 80]);
                    
                    % Set colormap and limits
                    colormap(AH(kk), utils_plot.redblue())
                    if bb <= 2
                        clim([-1 4]);
                        cbar_labels = {'-1', '1.5', '4'};
                        cticks = -1:2.5:4;
                    else
                        clim([1 2]);
                        cbar_labels = {'1', '1.5', '2'};
                        cticks = 1:0.5:2;
                    end
                    
                    % Add colorbar
                    h = utils_plot.createNarrowColorbar('vert');
                    % set(h,'ylim',clims);
                    set(h, 'YTick', cticks);
                    set(h, 'YTickLabels', cbar_labels);
                    colormap(h, utils_plot.redblue());
                    
                    title(sprintf('I0 in [%.2g, %.2g]', bins(bb), bins(bb+1)), ...
                          'FontSize', C.fonts.titleSize);
                    
                else
                    warning('No data points found for bin %d', bb);
                    scatter([0], [0], 5, 1, 'filled'); % Placeholder
                    axis([0 80 0 80]);
                end
                
            catch binError
                warning('Error processing bin %d: %s', bb, binError.message);
                axis([0 80 0 80]);
            end
            
            xlabel('r1 (Hz)', 'FontSize', C.fonts.labelSize);
            if kk == 3
                ylabel('r2 (Hz)', 'FontSize', C.fonts.labelSize);
            end
            
            pbaspect(C.figure.aspectRatio);
            kk = kk + 1;
        end
        
    catch ME
        warning('Error processing normalization patterns: %s', ME.message);
    end
    
    %% Theoretical panels
    n_range = [2.5 1 0.8];
    theta_range = [-0.32 -0.32];
    I0_range = [0.2 -0.4];
    
    % Panel showing different n values
    for nn = 1:length(n_range)
        iA = nn*3 + 1;
        if iA <= length(AH2)
            axes(AH2(iA));
            
            try
                n = n_range(nn);
                x = linspace(-0.5, 2, 101);
                y = x.^n .* (x > 0);
                plot(x, y, 'LineWidth', C.plot.lineWidth);
                title(sprintf('n=%.2g', n), 'FontSize', C.fonts.titleSize);
                pbaspect(C.figure.aspectRatio);
                
            catch ME
                warning('Error plotting theoretical panel %d: %s', nn, ME.message);
            end
        end
    end
    
    % Theoretical normalization patterns
    for kk = 1:2
        for nn = 1:length(n_range)
            iA = kk +1 + nn*3;
            if iA <= length(AH2)
                axes(AH2(iA));
                
                try
                    theta = theta_range(kk);
                    n = n_range(nn);
                    
                    [r1, r2] = meshgrid(0.1:0.1:90, 0.1:0.1:90);
                    I12 = (r1.^(1/n) + theta + r2.^(1/n) + theta)/1.5 + I0_range(kk);
                    r12 = (I12 - theta).^n .* (I12 >= theta);
                    NI = (r1 + r2) ./ r12;
                    NI(r12 <= 1) = NaN;
                    
                    imagesc(0.1:0.1:80, 0.1:0.1:80, NI);
                    axis xy;
                    axis([0 80 0 80]);
                    colormap(AH2(iA), utils_plot.redblue());
                    
                    C_range = max(abs(NI(:) - 1.5));
                    C_range = ceil(C_range * 10) / 10;
                    if C_range == 0
                        C_range = 0.1;
                    end
                    clim([1.5 - C_range, 1.5 + C_range]);
                    
                    h = utils_plot.createNarrowColorbar('vert');
                    set(h, 'YTick', [1 129 257]);
                    colormap(h, utils_plot.redblue());
                    
                    if nn == 1
                        title(sprintf('I0=%.2g', I0_range(kk)), 'FontSize', C.fonts.titleSize);
                    end
                    
                    xlabel('r1 (Hz)', 'FontSize', C.fonts.labelSize);
                    ylabel('r2 (Hz)', 'FontSize', C.fonts.labelSize);
                    pbaspect(C.figure.aspectRatio);
                    
                catch ME
                    warning('Error creating theoretical panel %d,%d: %s', kk, nn, ME.message);
                    imagesc(rand(80, 80)); % Placeholder
                    axis([0 80 0 80]);
                end
            end
        end
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_FFmodel', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_FFmodel:GeneralError', ...
          'Failed to generate feedforward model supplementary figure: %s', ME.message);
end

end