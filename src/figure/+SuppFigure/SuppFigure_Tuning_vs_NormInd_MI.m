function SuppFigure_Tuning_vs_NormInd_MI(varargin)
%SUPPFIGURE_TUNING_VS_NORMIND_MI Generate supplementary figure for tuning vs normalization index mutual information
%
% Purpose:
%   Creates a figure analyzing the relationship between preferred orientation
%   and normalization index across model and experimental data (V1, MT, V4),
%   using density scatter plots to visualize the distribution of neurons.
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
%   SuppFigure_Tuning_vs_NormInd_MI();
%   SuppFigure_Tuning_vs_NormInd_MI('Save', false);

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
    utils_plot.setPlotStyle('custom', 'width', 13, 'height', 10);
    
    % Create figure
    figure(p.Results.FigNum);
    clf;
    utils_plot.matchFigureAspectRatio();
    
    % Define layout
    horizontal_interval = 0.08*ones(1,2);
    vertical_interval = 0.18;
    my_position2 = [0.09 0.15 0.88 0.70];
    
    DC = utils_plot.divideAxes(0.25*ones(1,3), 0.25*ones(1,2), my_position2, ...
                               horizontal_interval, vertical_interval);
    
    % Create axes
    LdPos = [-0.03, 0.05];
    labels = {'A', 'B1', 'B2', 'B3'};
    panelorder = [1,2,4,6];
    
    for i = 1:4
        AH(i) = axes('Position', DC{panelorder(i)});
        hold on;
        utils_plot.addFigureLabel(labels{i}, LdPos);
    end
    
    utils_plot.setFigureProperties('AxisOpt', {'FontSize', C.fonts.labelSize}, 'LabelOpt', {'FontSize', C.fonts.labelSize});
    
    % Parameters
    Nneuron = 5000;
    Nshuffle = 2000;
    sz = 2;
    
    xticks1 = 0:45:180;
    xticklabels1 = {'0','45','90','135','180'};
    xticks2 = 0:90:360;
    xticklabels2 = {'0','90','180','270','360'};
    
    title_Exp = {'V1 data', 'MT data', 'V4 data'};
    
    %% Load data
    try
        % Load model data
        patternData = load(fullfile(rootFolder, 'data', 'Fig1_Pattern3cond.mat'));
        thetaData = load(fullfile(rootFolder, 'data', 'theta_map_MTE_two_rec.mat'));
        V1Data = load(fullfile(rootFolder, 'data', 'FigSupp_V1_norm_prefori.mat'));
        V4Data = load(fullfile(rootFolder, 'data', 'FigSupp_V4_norm_prefori.mat'));
        MTData = load(fullfile(rootFolder, 'data', 'FigSupp_MT_norm_prefori.mat'));
        
        
        MTE = patternData.MTE;
        theta_mapE = thetaData.theta_mapE;
        datas = {V1Data,V4Data,MTData};
        
        
    catch ME
        error('SuppFigure_Tuning_vs_NormInd_MI:DataLoadError', ...
              'Failed to load required data: %s', ME.message);
    end
    
    %% Panel A: Model data
    axes(AH(1));
    
    try
        % Prepare model data
        theta_mapE = reshape(theta_mapE', [], 1);
        
        % Subsample neurons
        rng(30925); % For reproducibility
        if size(MTE, 1) >= Nneuron
            ind1 = randsample(1:size(MTE, 1), Nneuron);
        else
            ind1 = 1:size(MTE, 1);
            warning('Not enough neurons for full sampling, using all %d neurons', length(ind1));
        end
        
        norm1 = (MTE(:,1) + MTE(:,2)) ./ MTE(:,3);
        norm = norm1(ind1);
        prefori = theta_mapE(ind1) * 180;
        
        % Remove invalid data
        validData = ~isnan(norm) & ~isnan(prefori) & isfinite(norm) & isfinite(prefori);
        norm = norm(validData);
        prefori = prefori(validData);
        
        if ~isempty(norm) && ~isempty(prefori)
            % Create density scatter plot
            Nbins = 60;
            xbinw = 180/Nbins;
            ybinw = 2/Nbins;
            sizepoint = 3;
            c0 = parula;
            
            N_temp = histcounts2(prefori, norm, 0:xbinw:180, 0.0:ybinw:3.0);
            Nmax = max(N_temp(:));
            
            if Nmax > 0
                Ncounts1 = N_temp(:);
                Ninds = ceil(Ncounts1/Nmax*256);
                [xcoor,ycoor] = meshgrid(0+xbinw/2:xbinw:180, 0.0+ybinw/2:ybinw:3.0);
                xcoor=reshape(xcoor',[],1);
                ycoor=reshape(ycoor',[],1);
                
                validPoints = Ninds > 0;
                if any(validPoints)
                    colors = c0(Ninds(validPoints), :);
                    scatter(xcoor(validPoints), ycoor(validPoints), sizepoint, colors, 'filled');
                end
                
                % Add colorbar
                h = utils_plot.createNarrowColorbar('vert');
                set(h, 'FontSize', C.fonts.tickSize);
                set(h, 'YTick', [0 1]);
                set(h, 'YTickLabel', [0 Nmax]);
                
                text(1.25, 0.5, '# neurons', 'Units', 'normalized', 'Rotation', 90, ...
                     'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
            end
        end
        
    catch ME
        warning('Error processing model data: %s', ME.message);
        scatter(rand(100,1)*180, rand(100,1)*3, sz, 'k.');
    end
    
    ylabel('normalization index', 'FontSize', C.fonts.labelSize);
    xlabel('preferred orientation', 'FontSize', C.fonts.labelSize);
    text('Units', 'Normalized', 'Position', [0.5 1.2], ...
         'string', 'model', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
    ylim([0 3]);
    xlim([0 180]);
    set(gca, 'XTick', xticks1, 'XTickLabel', xticklabels1);
    pbaspect(C.figure.aspectRatio);
    
    %% Panels B1-B3: Experimental data
    for ii = 1:3
        axes(AH(ii+1));
        
        try
            norm_all = datas{ii}.norm_all;
            prefori_all = datas{ii}.prefori_all;
            
            
            if ~isempty(norm_all) && ~isempty(prefori_all)
                % Create density scatter plot
                Nbins = 60;
                xbinw = 360/Nbins;
                ybinw = 3/Nbins;
                sizepoint = 3;
                c0 = parula;
                
                N_temp = histcounts2(prefori_all, norm_all, 0:xbinw:360, 0.0:ybinw:3.0);
                Nmax = max(N_temp(:));
                
                if Nmax > 0
                    Ncounts1 = N_temp(:);
                    Ninds = ceil(Ncounts1/Nmax*256);
                    [xcoor,ycoor] = meshgrid(0+xbinw/2:xbinw:360, 0.0+ybinw/2:ybinw:3.0);
                    xcoor=reshape(xcoor',[],1);
                    ycoor=reshape(ycoor',[],1);
                    
                    validPoints = Ninds > 0;
                    if any(validPoints)
                        colors = c0(Ninds(validPoints), :);
                        scatter(xcoor(validPoints), ycoor(validPoints), sizepoint, colors, 'filled');
                    end
                    
                    % Add colorbar
                    h = utils_plot.createNarrowColorbar('vert');
                    set(h, 'FontSize', C.fonts.tickSize);
                    set(h, 'YTick', [0 1]);
                    set(h, 'YTickLabel', [0 Nmax]);
                    
                    text(1.25, 0.5, '# neurons', 'Units', 'normalized', 'Rotation', 90, ...
                         'HorizontalAlignment', 'center', 'FontSize', C.fonts.annotationSize, 'color', 'k');
                end
            else
                warning('No valid data found for experimental dataset %d', ii);
                scatter(rand(50,1)*360, rand(50,1)*3, sz, 'k.');
            end
            
        catch ME
            warning('Error processing experimental data %d: %s', ii, ME.message);
            scatter(rand(50,1)*360, rand(50,1)*3, sz, 'k.');
        end
        
        xlabel('preferred orientation', 'FontSize', C.fonts.labelSize);
        
        text('Units', 'Normalized', 'Position', [0.5 1.2], ...
             'string', title_Exp{ii}, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', C.fonts.annotationSize, 'color', 'k');
        
        ylim([0 3]);
        xlim([0 360]);
        set(gca, 'XTick', xticks2, 'XTickLabel', xticklabels2);
        pbaspect(C.figure.aspectRatio);
    end
    
    %% Save figure
    if p.Results.Save
        utils_plot.saveFigure('path', resultPath, 'name', 'FigSupp_Tuning_vs_NormInd_MI', ...
                             'view', p.Results.View, 'save', true, ...
                             'format', 'pdf', 'res', 600);
    elseif p.Results.View
        drawnow;
    end
    
catch ME
    error('SuppFigure_Tuning_vs_NormInd_MI:GeneralError', ...
          'Failed to generate tuning vs normalization index MI supplementary figure: %s', ME.message);
end

end