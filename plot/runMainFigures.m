function runMainFigures(varargin)
% RUNALLFIGURES Generate all figures for the normalization paper
%
% Purpose:
%   Main script to generate all publication figures with proper error
%   handling and progress reporting.
%
% Inputs (optional):
%   'figures' - Array of figure numbers to generate (default: 1:6)
%   'verify'  - Logical, verify data files exist before running (default: true)
%   'parallel' - Logical, use parallel processing if available (default: false)
%
% Outputs:
%   All figures saved as PDFs in the results folder
%
% Usage Examples:
%   runMainFigures();                    % Generate all figures
%   runMainFigures('figures', [1, 3]);   % Generate only figures 1 and 3
%   runMainFigures('verify', false);     % Skip data verification

    % Parse inputs
    p = inputParser;
    addParameter(p, 'figures', 1:6, @(x) isnumeric(x) && all(x >= 1 & x <= 6));
    addParameter(p, 'verify', true, @islogical);
    addParameter(p, 'parallel', false, @islogical);
    parse(p, varargin{:});
    
    figureList = p.Results.figures;
    verifyData = p.Results.verify;
    useParallel = p.Results.parallel;
    
    % Initialize
    fprintf('\n=== FIGURE GENERATION STARTED ===\n');
    fprintf('Time: %s\n\n', datestr(now));
    
    % Setup paths
    C = figureConstants();
    thisFile = mfilename('fullpath');
    here     = fileparts(thisFile);          % .../+MainFigure
    rootFolder = fileparts(here);  % project root
    dataPath = fullfile(rootFolder, C.paths.dataFolder);
    resultPath = fullfile(rootFolder, C.paths.resultFolder);
    figurePath = fullfile(rootFolder, C.paths.figureFolder);
    utilsPath = fullfile(rootFolder, C.paths.utilsFolder);
    addpath(figurePath);
    addpath(utilsPath);
    
    % Create directories if needed
    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
        fprintf('Created results directory: %s\n', resultPath);
    end
    
    % Verify data files if requested
    if verifyData
        fprintf('Verifying data files...\n');
        [allPresent, missingFiles] = verifyDataFiles(dataPath);
        if ~allPresent
            fprintf('WARNING: Missing data files:\n');
            for i = 1:length(missingFiles)
                fprintf('  - %s\n', missingFiles{i});
            end
            response = input('Continue anyway? (y/n): ', 's');
            if ~strcmpi(response, 'y')
                fprintf('Aborted by user.\n');
                return;
            end
        else
            fprintf('All data files verified.\n\n');
        end
    end
    
    % Figure generation functions
    figureGenerators = {
        @MainFigure.generateFigure1,
        @MainFigure.generateFigure2,
        @MainFigure.generateFigure3,
        @MainFigure.generateFigure4,
        @MainFigure.generateFigure5,
        @MainFigure.generateFigure6
    };
    
    % Track timing and success
    nFigures = length(figureList);
    success = false(nFigures, 1);
    timings = zeros(nFigures, 1);
    
    % Generate figures
    if useParallel && license('test', 'Parallel_Computing_Toolbox')
        % Parallel execution
        fprintf('Using parallel processing...\n\n');
        parfor i = 1:nFigures
            figNum = figureList(i);
            [success(i), timings(i)] = generateSingleFigure(figNum, figureGenerators{figNum});
        end
    else
        % Sequential execution
        for i = 1:nFigures
            figNum = figureList(i);
            [success(i), timings(i)] = generateSingleFigure(figNum, figureGenerators{figNum});
        end
    end
    
    % Summary report
    fprintf('\n=== GENERATION SUMMARY ===\n');
    fprintf('Total time: %.2f seconds\n', sum(timings));
    fprintf('Successful: %d/%d figures\n', sum(success), nFigures);
    
    if sum(success) < nFigures
        fprintf('\nFailed figures:\n');
        failedIdx = find(~success);
        for i = 1:length(failedIdx)
            fprintf('  - Figure %d\n', figureList(failedIdx(i)));
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

function [allPresent, missingFiles] = verifyDataFiles(dataPath)
    % Verify all required data files exist
    
    % Complete list of required files
    requiredFiles = {
        % Figure 1
        'theta_map_MTE_two_rec.mat'
        'weight_recurrent.mat'
        'Fig1_GaborImageDemonstration.mat'
        'Fig1_Pattern3cond.mat'
        'Fig1_NormIndDistribution.mat'
        'Fig_ExpData_V4_NormIndDistribution.mat'
        'Fig_ExpData_MT_NormIndDistribution.mat'
        % Figure 2
        'Fig2_Norm1Norm2Corr.mat'
        'Fig2_DiagonalCorrDecreaseWithAveNMI.mat'
        'Fig2_AntiDiagonalCorrDecreaseWithDiffNMI.mat'
        'Fig2_TuningSimilarityControl.mat'
        'Fig_ExpData_MT_Norm1Norm2Corr.mat'
        'Fig_ExpData_MT_DiagonalCorrDecreaseWithAveNMI_Width_0d2.mat'
        'Fig_ExpData_MT_AntiDiagonalCorrDecreaseWithDiffNMI.mat'
        'Fig_ExpData_MT_TuningSimilarityControl_weighted.mat'
        'Fig_ExpData_V4_Norm1Norm2Corr.mat'
        'Fig_ExpData_V4_DiagonalCorrDecreaseWithAveNMI_Width_0d2.mat'
        'Fig_ExpData_V4_AntiDiagonalCorrDecreaseWithDiffNMI.mat'
        'Fig_ExpData_V4_TuningSimilarityControl_weighted.mat'
        % Figure 3
        'Fig3_NormIndFiringRateCurrent.mat'
        'Fig3_CurrentCovarianceFixOneAxis.mat'
        'Fig3_CurrentCovarianceDiagonal.mat'
        % Figure 4
        'Fig4_fr_con1_con2.mat'
        % Figure 5
        'Fig5_FI_Con_Default_RandomNum_Arithmetic.mat'
        'Fig5_FI_Ori_Default_RandomNum_Arithmetic.mat'
        % Figure 6
        'Fig6_FRpattern_InDegree.mat'
        'Fig6_FI_Con.mat'
        'Fig6_FI_Ori.mat'
        'Fig6_manifold_Default.mat'
        'Fig6_manifold_InDegree.mat'
        % Shared colormap
        'mycolormap2.mat'
    };
    
    % Check each file
    missingFiles = {};
    for i = 1:length(requiredFiles)
        filepath = fullfile(dataPath, requiredFiles{i});
        if ~exist(filepath, 'file')
            missingFiles{end+1} = requiredFiles{i};
        end
    end
    
    allPresent = isempty(missingFiles);
end