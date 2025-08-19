function [Wrr, Wrf, degree_params] = load_indegree_connectivity(options, params)
% LOAD_INDEGREE_CONNECTIVITY - Load connectivity with degree information
%
% H1 Line: Load connectivity matrices with heterogeneous degrees
%
% Syntax:
%   [Wrr, Wrf, degree_params] = load_indegree_connectivity(options, params)
%
% Description:
%   Loads connectivity matrices along with degree distribution
%   information for heterogeneous network simulations.
%
% Inputs:
%   options - Simulation options with W_fname
%   params  - Network parameters
%
% Outputs:
%   Wrr           - Recurrent connectivity
%   Wrf           - Feedforward connectivity
%   degree_params - Structure with degree information
%
% Example:
%   [Wrr, Wrf, deg_p] = load_indegree_connectivity(options, params);

    if ~exist(options.W_fname, 'file')
        error('SpikingSim:FileNotFound', ...
              'Weight file not found: %s', options.W_fname);
    end
    
    % Load all connectivity data
    data = load(options.W_fname, 'Wrr', 'Wrf', ...
                'outdegreeEE', 'outdegreeIE', 'outdegreeEI', 'outdegreeII', ...
                'outdegreeEX', 'outdegreeIX', 'seqindX', 'seqindE', 'seqindI');
    
    Wrr = data.Wrr;
    Wrf = data.Wrf;
    
    % Package degree information
    degree_params = struct();
    degree_params.outdegreeEE = int32(data.outdegreeEE);
    degree_params.outdegreeEI = int32(data.outdegreeEI);
    degree_params.outdegreeIE = int32(data.outdegreeIE);
    degree_params.outdegreeII = int32(data.outdegreeII);
    degree_params.outdegreeEX = int32(data.outdegreeEX);
    degree_params.outdegreeIX = int32(data.outdegreeIX);
    degree_params.seqindX = int32(data.seqindX);
    degree_params.seqindE = int32(data.seqindE);
    degree_params.seqindI = int32(data.seqindI);
    
    % Calculate maximum degrees
    degree_params.Kemax = max(data.outdegreeEE + data.outdegreeIE);
    degree_params.Kimax = max(data.outdegreeEI + data.outdegreeII);
    degree_params.Kxmax = max(data.outdegreeEX + data.outdegreeIX);
    
    fprintf('Loaded in-degree matched connectivity from: %s\n', options.W_fname);
    fprintf('Max degrees - E: %d, I: %d, X: %d\n', ...
            degree_params.Kemax, degree_params.Kimax, degree_params.Kxmax);
end
