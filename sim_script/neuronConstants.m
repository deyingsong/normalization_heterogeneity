function param = neuronConstants()
%NEURONCONSTANTS  One-stop constants for the two-layer visual cortex model.
% Syntax
%   P = neuronConstants();
%
% Description
%   Returns a struct P with fields describing geometry, neuron counts,
%   V1 (input) statistics, and V4/MT (recurrent spiking) connectivity and
%   synaptic parameters, consolidated from the provided manuscript pages.
%   Use this once at startup and pass P around instead of scattering
%   numbers through the code.
%
% Output
%   P : struct with subfields
%       .geom, .V1, .V4MT, .conn, .syn
%
%
% -------------------------------------------------------------------------

P = struct();

%% Geometry / domain
P.geom.domain = [0 2; 0 1];       % Γ = [x_min x_max; y_min y_max]

%% Population sizes
P.Nx = 5e3;                % V1
P.Ne = 2e4;                % V4/MT excitatory
P.Ni = 5e3;                % V4/MT inhibitory


%% V1 model neuron (receptive fields & drive)
P.V1.rf.sigma   = 0.2;            % Gaussian envelope std (Gabor)
P.V1.rf.lambda  = 0.6;            % wavelength of sinusoid
P.V1.rf.centers = [0.5 0.5;       % center for V11
                   1.5 0.5];      % center for V12 (both in Γ coordinates)

% Orientation map 
P.V1.map.nColumns = 30;
P.V1.map.Lambda   = 0.125;        % average column spacing
P.V1.map.ljValues = [+1, -1];     % random ±1 per column
P.V1.map.phiRange = [0, 2*pi];    % φ_j ~ U[0,2π]

% Poisson drive statistics
P.V1.rates.mean_when_image = 10;  % Hz, average rate in the presence of Gabor image
P.V1.rates.spontaneous     = 5;   % Hz, baseline when image absent

%% Connection geometry (distance kernels) and selection rules
% Wrapped Gaussian width (periodic on Γ)
P.conn.sigma.recurrent_exc = 0.2; % σ_e for V4/MT recurrent
P.conn.sigma.recurrent_inh = 0.2; % σ_i for V4/MT recurrent
P.conn.sigma.feedforward   = 0.1; % σ_ffwd from V1 → V4/MT

% Fractions and tuning rule for excitatory projections
P.conn.excit.distance_fraction = 0.85;     % 85% distance-dependent
P.conn.excit.tuned_fraction    = 0.15;     % 15% similarly-tuned (no spatial dependence)
P.conn.excit.tuned_cos_thresh  = 0.6;      % cos(θ_i − θ_j) ≥ 0.6

% Inhibitory projections ignore tuning
P.conn.inhib.distance_only = true;

% Helper: wrapped Gaussian (periodized separable 2D). Use finite sums in practice.
P.conn.wrapGauss = @(dx,dy,sigma) (1/(2*pi*sigma^2)) * ...
    (sum(exp(-((dx + ( -3:3 )*2).^2) /(2*sigma^2))) .* ...
     sum(exp(-((dy + ( -3:3 )*1).^2) /(2*sigma^2)))); 
% (Implementation note: the above uses a small periodization window; adjust if needed.)

%% Synaptic strengths, probabilities and out-degrees (V4/MT and feedforward)
% Recurrent synaptic weights within V4/MT (mV), scaled by 1/sqrt(N)
P.syn.scale_with_inv_sqrtN = true;
P.syn.J.ee =  +80;    % mV
P.syn.J.ei = -240;    % mV
P.syn.J.ie =  +40;    % mV
P.syn.J.ii = -300;    % mV

% Mean connection probabilities (recurrent, V4/MT)
P.conn.pbar.ee = 0.01;
P.conn.pbar.ei = 0.04;
P.conn.pbar.ie = 0.03;
P.conn.pbar.ii = 0.04;

% Out-degrees (recurrent, V4/MT)
P.conn.Kout.ee = 200;
P.conn.Kout.ei = 800;
P.conn.Kout.ie = 150;
P.conn.Kout.ii = 200;

% Feedforward V1 → V4/MT synapses
P.syn.J.eF = 160;     % mV, onto E cells
P.syn.J.iF = 140;     % mV, onto I cells
P.conn.pbar.eF = 0.05;
P.conn.pbar.iF = 0.05;
P.conn.Kout.eF = 1000;
P.conn.Kout.iF =  250;

% Multiple synapses from one pre to one post are allowed
P.conn.allow_multisyn = true;

%% Package and return
param = P;
end
