function params = initialize_network_parameters()
% INITIALIZE_NETWORK_PARAMETERS - Set up default network parameters
%
% H1 Line: Initialize default parameters for spiking neural network
%
% Syntax:
%   params = initialize_network_parameters()
%
% Description:
%   Creates a structure containing all default parameters for the
%   spiking neural network simulation including network size, 
%   connectivity, neuron properties, and simulation settings.
%
% Outputs:
%   params - Structure containing all network parameters
%
% Example:
%   params = initialize_network_parameters();

    % Network dimensions
    params.dim = '2D';
    params.Nx1 = 100;  % Feedforward layer dimension
    params.Ne1 = 200;  % Excitatory layer dimension  
    params.Ni1 = 100;  % Inhibitory layer dimension
    
    % Calculate total populations
    params.Ne = params.Ne1 * params.Ne1 / 2;  % Total excitatory neurons
    params.Ni = params.Ni1 * params.Ni1 / 2;  % Total inhibitory neurons
    params.Nx = params.Nx1 * params.Nx1 / 2;  % Total feedforward neurons
    params.N = params.Ne + params.Ni;         % Total recurrent neurons
    
    % Simulation parameters
    params.T = 20000;        % Total simulation time (ms)
    params.dt = 0.05;        % Time step (ms)
    params.Tburn = 10;       % Burn-in period (ms)
    params.maxrate = 0.05;   % Maximum average firing rate (kHz)
    
    % Static currents
    params.inE = 0;  % External current to excitatory neurons
    params.inI = 0;  % External current to inhibitory neurons
    params.Iapp = [params.inE; params.inI];
    
    % Connection widths (spatial extent of connections)
    params.sigmaRX = 0.1 * ones(2, 1);  % Feedforward to recurrent
    params.sigmaRR = 0.2 * ones(2, 2);  % Recurrent to recurrent
    
    % Connection probabilities
    params.Prr = [0.01, 0.04;  % E->E, E->I
                  0.03, 0.04]; % I->E, I->I
    params.Prx = [0.05; 0.05]; % X->E, X->I
    
    % Connection strengths (before scaling)
    params.Jr = [80, -240;   % E->E, E->I
                 40, -300];  % I->E, I->I
    params.Jx = [160, 140, 160, 140];  % Feedforward strengths
    
    % Synaptic time constants
    params.taudsyn = [5, 100;   % Decay time constants (X)
                      5, 100;   % Decay time constants (E)
                      8, 100];  % Decay time constants (I)
    params.taursyn = [1, 2;     % Rise time constants (X)
                      1, 2;     % Rise time constants (E)
                      1, 2];    % Rise time constants (I)
    params.Psyn = [0.2, 0.8;    % Synapse type percentages (X)
                   1.0, 0.0;    % Synapse type percentages (E)
                   1.0, 0.0];   % Synapse type percentages (I)
    
    % EIF neuron parameters
    params.gl = [1/15, 1/10];        % Leak conductance (E, I)
    params.Cm = [1, 1];              % Membrane capacitance (E, I)
    params.vlb = [-100, -100];       % Lower bound voltage (E, I)
    params.vth = [-10, -10];         % Threshold voltage (E, I)
    params.DeltaT = [2, 0.5];        % Slope factor (E, I)
    params.vT = [-50, -50];          % Spike threshold (E, I)
    params.vre = [-65, -65];         % Reset voltage (E, I)
    params.tref = [1.5, 0.5];        % Refractory period (E, I)
    
    % Recording settings
    params.nrecordE0 = 1;  % Number of E neurons to record
    params.nrecordI0 = 1;  % Number of I neurons to record
    
    % Additional parameters for specific simulations
    params.sigma_current = 6.8;  % Current noise standard deviation
    params.sigma_n = 3.5;        % Input noise parameter
    params.tau_n = 40;           % Input noise time constant
    params.NI = 450;             % Number of input channels
end
