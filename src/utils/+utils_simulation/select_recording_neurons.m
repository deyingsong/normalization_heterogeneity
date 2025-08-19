function Irecord = select_recording_neurons(params)
% SELECT_RECORDING_NEURONS - Choose neurons to record
%
% H1 Line: Select random neurons for detailed recording
%
% Syntax:
%   Irecord = select_recording_neurons(params)
%
% Description:
%   Randomly selects excitatory and inhibitory neurons for
%   recording synaptic currents and membrane potentials.
%
% Inputs:
%   params - Network parameters with nrecordE0 and nrecordI0
%
% Outputs:
%   Irecord - Vector of neuron indices to record
%
% Example:
%   Irecord = select_recording_neurons(params);

    exc_neurons = randi(params.Ne, 1, params.nrecordE0);
    inh_neurons = randi(params.Ni, 1, params.nrecordI0) + params.Ne;
    Irecord = [exc_neurons, inh_neurons];
end
