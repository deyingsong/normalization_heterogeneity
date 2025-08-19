function filters = load_filters(filename)
% LOAD_FILTERS - Load visual filter matrix
%
% H1 Line: Load filter matrix from MAT file
%
% Syntax:
%   filters = load_filters(filename)
%
% Description:
%   Loads the visual filter matrix used for input generation.
%   Handles missing files gracefully.
%
% Inputs:
%   filename - MAT file containing filter matrix
%
% Outputs:
%   filters - Structure containing filter matrix
%
% Example:
%   filters = load_filters('V1filterRecSig0d2Lam0d6.mat');

    if ~exist(filename, 'file')
        error('SpikingSim:FilterFileNotFound', ...
              'Filter file not found: %s', filename);
    end
    
    filters = load(filename);
    
    if ~isfield(filters, 'F')
        error('SpikingSim:MissingFilterMatrix', ...
              'Filter file does not contain matrix F');
    end
end
