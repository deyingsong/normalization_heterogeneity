function results = analyze_simulation_results(s2, params, elapsed_time)
% ANALYZE_SIMULATION_RESULTS - Calculate simulation statistics
%
% H1 Line: Analyze spike data and compute network statistics
%
% Syntax:
%   results = analyze_simulation_results(s2, params, elapsed_time)
%
% Description:
%   Computes firing rates, synchrony measures, and other statistics
%   from the simulation output.
%
% Inputs:
%   s2           - Spike train output [time; neuron_id]
%   params       - Network parameters
%   elapsed_time - Simulation runtime
%
% Outputs:
%   results - Structure containing analysis results
%
% Example:
%   results = analyze_simulation_results(s2, params, 10.5);

    % Find valid spikes
    if isempty(s2)
        warning('SpikingSim:NoSpikes', 'No spikes generated');
        results = struct('nuSim', [0, 0], 'elapsed_time', elapsed_time);
        return;
    end
    
    End = find(s2(2,:) == 0, 1) - 1;
    if isempty(End)
        End = size(s2, 2);
    end
    
    % Calculate firing rates after burn-in
    exc_spikes = s2(1,1:End) > params.Tburn & s2(2,1:End) <= params.Ne;
    inh_spikes = s2(1,1:End) > params.Tburn & s2(2,1:End) > params.Ne;
    
    T_analysis = params.T - params.Tburn;
    
    results.nuSim(1) = 1000 * nnz(exc_spikes) / (params.Ne * T_analysis); % Hz
    results.nuSim(2) = 1000 * nnz(inh_spikes) / (params.Ni * T_analysis); % Hz
    
    % Print results
    fprintf('\n=== Simulation Results ===\n');
    fprintf('Average firing rates:\n');
    fprintf('  Excitatory: %.2f Hz\n', results.nuSim(1));
    fprintf('  Inhibitory: %.2f Hz\n', results.nuSim(2));
    fprintf('Elapsed time: %.2f seconds\n', elapsed_time);
    fprintf('========================\n\n');
    
    % Additional statistics
    results.elapsed_time = elapsed_time;
    results.total_spikes = End;
    results.spikes_per_neuron = End / params.N;
    
    % Check for pathological states
    if results.nuSim(1) > 100
        warning('SpikingSim:HighRate', ...
                'Excitatory rate unusually high (%.1f Hz)', results.nuSim(1));
    end
    if results.nuSim(1) < 0.1
        warning('SpikingSim:LowRate', ...
                'Excitatory rate very low (%.3f Hz)', results.nuSim(1));
    end
end
