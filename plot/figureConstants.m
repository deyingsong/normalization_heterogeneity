function C = figureConstants()
% FIGURECONSTANTS Define global constants for all figures
%
% Purpose:
%   Centralizes all constants used across figures to maintain consistency
%   and avoid magic numbers in the code.
%
% Outputs:
%   C - Structure containing all constants organized by category
%
% Usage Example:
%   C = figureConstants();
%   colors = C.colors.excitatory;

    % Network parameters
    C.network.Ne = 2e4;           % Number of excitatory neurons
    C.network.Ni = 5e3;           % Number of inhibitory neurons
    C.network.N = C.network.Ne + C.network.Ni;  % Total neurons
    C.network.Tburn = 1000;       % Burn-in time
    
    % Color schemes
    C.colors.excitatory = [202, 0, 32] / 255;      % Red for excitatory
    C.colors.inhibitory = [5, 113, 176] / 255;     % Blue for inhibitory
    C.colors.feedforward = [244, 165, 130] / 255;  % Light red
    C.colors.gray = [0.7, 0.7, 0.7];               % Standard gray
    
    % Figure 3 specific colors
    C.colors.fig3_palette = [
        202, 0, 32;       % Deep red
        228, 59, 126;     % Pink
        247, 247, 247;    % White
        0, 198, 255;      % Light blue
        5, 113, 176       % Deep blue
    ] / 255;
    
    % Figure 4 copper colors
    C.colors.copper_palette = copper(6);
    
    % Hot colormap for Figure 5
    C.colors.hot_palette = hot(8);
    
    % Plot parameters
    C.plot.markerSize = 4;
    C.plot.lineWidth = 1;
    C.plot.markerLineWidth = 1.5;
    C.plot.scatterPointSize = 3;
    C.plot.alphaValue = 0.2;  % Transparency for scatter plots
    
    % Figure dimensions
    C.figure.width = 14;      % Default figure width
    C.figure.height = 10;     % Default figure height
    C.figure.aspectRatio = [1, 1, 1];  % Default aspect ratio
    
    % Axis limits and bins
    C.analysis.normIndexBins = 0:0.1:3;
    C.analysis.normIndexRange = [0, 3];
    C.analysis.firingRateRange = [0, 90];
    C.analysis.contrastRange = [0, 1];
    
    % Statistical parameters
    C.stats.nBins = 60;
    C.stats.nSamples = 20;
    C.stats.binWidth = 0.05;
    
    % File paths (relative to project root)
    C.paths.dataFolder = 'data';
    C.paths.resultFolder = 'results';
    C.paths.utilsFolder = 'src/utils';
    C.paths.figureFolder = 'src/figure';
    
    % Font sizes
    C.fonts.titleSize = 9;
    C.fonts.labelSize = 8;
    C.fonts.annotationSize = 8;
    C.fonts.tickSize = 8;
    
    % Legend and label positions
    C.positions.legendOffset = [-0.04, 0.04];
    C.positions.titleHeight = 1.25;
    C.positions.subtitleHeight = 1.12;
end
