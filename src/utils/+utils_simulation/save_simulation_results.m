function save_simulation_results(filename, s2, results, params, options)
% SAVE_SIMULATION_RESULTS - Save simulation output to file
%
% H1 Line: Save simulation data and results to MAT file
%
% Syntax:
%   save_simulation_results(filename, s2, results, params, options)
%
% Description:
%   Saves spike data, analysis results, and parameters to a MAT file
%   with optional compression for large datasets.
%
% Inputs:
%   filename - Output filename
%   s2      - Spike data
%   results - Analysis results
%   params  - Network parameters
%   options - Simulation options
%
% Example:
%   save_simulation_results('output.mat', s2, results, params, options);

    fprintf('Saving results to: %s\n', filename);
    
    % Remove zero entries from spike data if requested
    if options.saveS2
        s2 = s2(:, s2(2,:) ~= 0);
    end
    
    % Create output structure
    output = struct();
    output.T = params.T;
    output.nuSim = results.nuSim;
    
    % Save basic results
    save(filename, '-struct', 'output', '-v7.3');
    
    % Optionally save spike data
    if options.saveS2
        save(filename, 's2', '-append');
    end
    
    % Optionally save parameters
    if options.saveParam
        save(filename, 'params', '-append');
    end
    
    % Optionally save population rates
    if options.saveRm
        re2 = hist(s2(1, s2(2,:) <= params.Ne & s2(2,:) > 0), 1:params.T) ...
              / params.Ne * 1e3;
        ri2 = hist(s2(1, s2(2,:) > params.Ne), 1:params.T) ...
              / params.Ni * 1e3;
        save(filename, 're2', 'ri2', '-append');
    end
    
    fprintf('Results saved successfully\n');
end
