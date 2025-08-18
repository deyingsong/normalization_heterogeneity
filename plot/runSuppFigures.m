function runSuppFigures(varargin)
% RUNALLFIGURES Generate all figures for the normalization paper
%
% Purpose:
%   Main script to generate all publication figures with proper error
%   handling and progress reporting.
%
% Inputs (optional):
%   'figures' - Array of figure numbers to generate (default: 1:6)
%   'parallel' - Logical, use parallel processing if available (default: false)
%
% Outputs:
%   All figures saved as PDFs in the results folder
%
% Usage Examples:
%   runSuppFigures();                    % Generate all figures

    % Parse inputs
    p = inputParser;

    % Define the allowed figure names
    validFigures = { ...
        'BasicResponse', 'Broad_InDegree', 'BroadWeight', ...
        'Current_Connection_NormInd', 'DynamicalRegime_GaussianNoise', ...
        'DynamicalRegime_wave', 'FFmodel_I0', 'FFmodel', ...
        'fI_curve_control', 'I_activity', 'PrefOri_StimOri', ...
        'RandomNet', 'SSN', 'Superimposed', 'V1', ...
        'VarOriPair', 'WeakCoupling_EIbalance'};    
    % Add parameter with default = all figure names
    addParameter(p, 'figures', validFigures, ...
    @(x) all(ismember(x, validFigures)));
    addParameter(p, 'parallel', false, @islogical);

    % Parse inputs
    
    parse(p, varargin{:});    
    figureList = p.Results.figures;
    useParallel = p.Results.parallel;
    
    % Initialize
    fprintf('\n=== FIGURE GENERATION STARTED ===\n');
    fprintf('Time: %s\n\n', datestr(now));
    
    % Setup paths
    C = figureConstants();
    dataPath = fullfile(pwd, C.paths.dataFolder);
    resultPath = fullfile(pwd, C.paths.resultFolder);
    figurePath = fullfile(pwd, C.paths.figureFolder);
    utilsPath = fullfile(pwd, C.paths.utilsFolder);
    addpath(figurePath);
    addpath(utilsPath);
    
    % Create directories if needed
    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
        fprintf('Created results directory: %s\n', resultPath);
    end
    
    
    % Figure generation functions
    % figureGenerators = {
    %     @SuppFigure.SuppFigure_V1,
    % };
    figureGenerators = containers.Map( ...
        validFigures, ...
        cellfun(@(n) str2func(['SuppFigure.SuppFigure_' n]), validFigures, 'UniformOutput', false));
    
    if isstring(figureList); figureList = cellstr(figureList); end
    if ischar(figureList);   figureList = {figureList};        end

    % Track timing and success
    nFigures = length(figureList);
    success = false(nFigures, 1);
    timings = zeros(nFigures, 1);
    
    % Generate figures
    if useParallel && license('test', 'Parallel_Computing_Toolbox')
        % Parallel execution
        fprintf('Using parallel processing...\n\n');
        parfor i = 1:nFigures
            name = figureList{i};
            if ~isKey(figureGenerators, name)
                warning('No generator found for "%s". Skipping.', name);
                continue
            end
            fh = figureGenerators(name);  % function handle, e.g., @SuppFigure_V1
        
            t0 = tic;
            try
                fh();                     % run the figure function with no args
                success(i) = true;
            catch ME
                warning('Figure "%s" failed: %s', name, ME.message);
                success(i) = false;
            end
            timings(i) = toc(t0);
        end
    else
        % Sequential execution
        for i = 1:nFigures
            name = figureList{i};
            if ~isKey(figureGenerators, name)
                warning('No generator found for "%s". Skipping.', name);
                continue
            end
            fh = figureGenerators(name);  % function handle, e.g., @SuppFigure_V1
        
            t0 = tic;
            try
                fh();                     % run the figure function with no args
                success(i) = true;
            catch ME
                warning('Figure "%s" failed: %s', name, ME.message);
                success(i) = false;
            end
            timings(i) = toc(t0);
        end
    end
    
    % Summary report
    fprintf('\n=== GENERATION SUMMARY ===\n');
    fprintf('Total time: %.2f seconds\n', sum(timings));
    fprintf('Successful: %d/%d figures\n', sum(success), nFigures);
    
    if sum(success) < nFigures
        fprintf('\nFailed figures:\n');
        failedIdx = find(~success);
        for i = 1:numel(failedIdx)
            fprintf('  - %s\n', figureList{failedIdx(i)});
        end
    end
    
    fprintf('\nResults saved in: %s\n', resultPath);
    fprintf('=== COMPLETE ===\n\n');
end

function [success, timing] = generateSingleFigure(figNum, generatorFunc)
    % Generate a single figure with error handling
    
    fprintf('Generating Figure %d... ', figNum);
    startTime = tic;
    success = false;
    
    try
        generatorFunc();
        timing = toc(startTime);
        success = true;
        fprintf('Done (%.2f s)\n', timing);
    catch ME
        timing = toc(startTime);
        fprintf('FAILED\n');
        fprintf('  Error: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  Location: %s (line %d)\n', ...
                ME.stack(1).file, ME.stack(1).line);
        end
    end
end

