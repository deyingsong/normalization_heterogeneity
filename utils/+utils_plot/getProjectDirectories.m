function directories = getProjectDirectories(varargin)
% GETPROJECTDIRECTORIES Get or set project directory paths
%
% Syntax:
%   dirs = getProjectDirectories()
%   dirs = getProjectDirectories('Name', 'path', ...)
%
% Description:
%   Manages project-specific directory paths with platform awareness.
%
% Inputs:
%   Name-value pairs of directory names and paths
%
% Outputs:
%   directories - Structure containing all directory paths
%
% Example:
%   dirs = getProjectDirectories('Data', '/path/to/data', 'Results', '/path/to/results');

    % Global storage
    persistent projectDirs;
    
    % Initialize if needed
    if isempty(projectDirs)
        projectDirs = initializeDefaultDirectories();
    end
    
    % Parse new directories if provided
    if nargin > 0
        newDirs = parseNameValuePairs(varargin{:});
        fieldNames = fieldnames(newDirs);
        
        for i = 1:length(fieldNames)
            % Validate path
            path = newDirs.(fieldNames{i});
            if ~ischar(path) && ~isstring(path)
                warning('getProjectDirectories:InvalidPath', ...
                        'Path for %s must be a string', fieldNames{i});
                continue;
            end
            
            % Normalize path
            path = char(path);
            if ~endsWith(path, filesep)
                path = [path, filesep];
            end
            
            % Store
            projectDirs.(fieldNames{i}) = path;
            
            % Check existence
            if ~exist(path, 'dir')
                warning('getProjectDirectories:PathNotFound', ...
                        'Directory does not exist: %s', path);
            end
        end
    end
    
    directories = projectDirs;
    
    % --- Helper Function ---
    function dirs = initializeDefaultDirectories()
        % Set platform-specific defaults
        if ispc
            rootDir = 'C:\Projects\';
        elseif ismac
            rootDir = '~/Projects/';
        else
            rootDir = '~/projects/';
        end
        
        % Expand user directory
        if startsWith(rootDir, '~')
            rootDir = fullfile(getenv('HOME'), rootDir(2:end));
        end
        
        % Default structure
        dirs = struct();
        dirs.Root = rootDir;
        dirs.Data = fullfile(rootDir, 'Data', filesep);
        dirs.Results = fullfile(rootDir, 'Results', filesep);
        dirs.Figures = fullfile(rootDir, 'Figures', filesep);
        dirs.Code = fullfile(rootDir, 'Code', filesep);
    end
end
