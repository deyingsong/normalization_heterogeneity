function V0 = initialize_membrane_potentials(params)
% INITIALIZE_MEMBRANE_POTENTIALS - Set initial voltages
%
% H1 Line: Initialize membrane potentials randomly
%
% Syntax:
%   V0 = initialize_membrane_potentials(params)
%
% Description:
%   Generates random initial membrane potentials between reset
%   and threshold values.
%
% Inputs:
%   params - Network parameters
%
% Outputs:
%   V0 - Initial voltage vector
%
% Example:
%   V0 = initialize_membrane_potentials(params);

    V0 = (params.V0max - params.V0min) .* rand(params.N, 1) + params.V0min;
end
