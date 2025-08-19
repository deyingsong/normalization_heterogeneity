function sx = gen_fisher_info_spikes(fr, T, F1, F2, NI, sigma_n, tau_n, ConPerc)
% GENERATE_FISHER_INFO_SPIKES - Generate spikes for Fisher information analysis
%
% H1 Line: Create spike trains with ON/OFF structure and contrast modulation
%
% Syntax:
%   sx = generate_fisher_info_spikes(fr, T, F1, F2, NI, sigma_n, tau_n, ConPerc)
%
% Description:
%   Generates spike trains specifically designed for Fisher information
%   analysis. Implements ON/OFF cell populations with contrast modulation
%   during specific time windows. The stimulus structure includes:
%   - Baseline period (0-1 ms)
%   - Stimulus period with contrast (1-200 ms)
%   - Post-stimulus period (201-500 ms)
%   This pattern repeats throughout the simulation.
%
% Inputs:
%   fr      - Baseline firing rates for all neurons (Hz)
%   T       - Total simulation time (ms)
%   F1      - Filter matrix for ON cells (2500 x 225)
%   F2      - Filter matrix for OFF cells (2500 x 225)
%   NI      - Number of noise inputs
%   sigma_n - Noise standard deviation
%   tau_n   - Noise correlation time constant (ms)
%   ConPerc - Contrast percentage (0-1)
%
% Outputs:
%   sx - Spike train matrix (2 x nspikes)
%        sx(1,:) - spike times
%        sx(2,:) - neuron indices
%
% Example:
%   sx = generate_fisher_info_spikes(fr, 20000, F1, F2, 450, 3.5, 40, 0.5);
%
% See also: generate_input_spikes

    % Validate inputs
    validate_fisher_parameters(fr, T, F1, F2, NI, ConPerc);
    
    Nx = length(fr);
    dt = 1;  % 1 ms time step
    
    % Pre-allocate spike array
    expected_spikes = sum(fr) * T / 1000 * 1.5;
    sx = zeros(2, round(expected_spikes * 2));
    nspks = 0;
    
    % Initialize noise processes
    noise1 = sigma_n / sqrt(2 * tau_n) * randn(NI/2, 1);
    noise2 = sigma_n / sqrt(2 * tau_n) * randn(NI/2, 1);
    
    % Define ON and OFF cell indices
    [ind1, ind2] = create_on_off_indices();
    
    % Initialize firing rate vector
    FR = zeros(5000, 1);
    
    % Check dimensions
    if length(ind1) ~= size(F1, 1) || length(ind2) ~= size(F2, 1)
        error('FisherSpikes:DimensionMismatch', ...
              'Filter matrices dimensions do not match population indices');
    end
    
    % Main simulation loop - process in 500ms blocks
    n_blocks = T / 500;
    if floor(n_blocks) ~= n_blocks
        warning('FisherSpikes:TimeNotDivisible', ...
                'T is not divisible by 500, truncating to %d ms', ...
                floor(n_blocks) * 500);
        n_blocks = floor(n_blocks);
    end
    
    for block = 1:n_blocks
        sx = process_fisher_block(block, sx, nspks, fr, F1, F2, ...
                                 noise1, noise2, ind1, ind2, ...
                                 FR, ConPerc, sigma_n, tau_n, dt);
        nspks = nnz(sx(1, :) > 0);
    end
    
    % Trim and sort final spike train
    sx = sx(:, sx(1, :) > 0 & sx(1, :) <= T);
    [~, J] = sort(sx(1, :));
    sx = sx(:, J);
end

function validate_fisher_parameters(fr, T, F1, F2, NI, ConPerc)
% VALIDATE_FISHER_PARAMETERS - Validate Fisher info spike parameters
%
% H1 Line: Check parameters for Fisher information spike generation

    % Basic parameter validation
    validate_spike_parameters(fr, T);
    
    % Check filter matrices
    if size(F1, 1) ~= 2500 || size(F1, 2) ~= 225
        error('FisherSpikes:InvalidF1', ...
              'F1 must be 2500x225 (got %dx%d)', size(F1, 1), size(F1, 2));
    end
    
    if size(F2, 1) ~= 2500 || size(F2, 2) ~= 225
        error('FisherSpikes:InvalidF2', ...
              'F2 must be 2500x225 (got %dx%d)', size(F2, 1), size(F2, 2));
    end
    
    % Check NI
    if NI ~= 450
        warning('FisherSpikes:UnexpectedNI', ...
                'NI=%d, expected 450 for standard configuration', NI);
    end
    
    % Check contrast
    if ConPerc < 0 || ConPerc > 1
        error('FisherSpikes:InvalidContrast', ...
              'ConPerc must be between 0 and 1 (got %.2f)', ConPerc);
    end
    
    % Check that T is reasonable for Fisher info
    if T < 500
        error('FisherSpikes:ShortDuration', ...
              'T must be at least 500ms for one stimulus block');
    end
end

function [ind1, ind2] = create_on_off_indices()
% CREATE_ON_OFF_INDICES - Generate indices for ON and OFF populations
%
% H1 Line: Create index vectors for ON and OFF cell populations
%
% Syntax:
%   [ind1, ind2] = create_on_off_indices()
%
% Description:
%   Creates index vectors that map the 50x50 ON and OFF populations
%   to the full 100x100 grid layout.
%
% Outputs:
%   ind1 - Indices for ON cells (2500x1)
%   ind2 - Indices for OFF cells (2500x1)

    % ON cells: first 50 columns of each row
    ind1 = zeros(50*50, 1);
    for i = 1:50
        for j = 1:50
            ind1((i-1)*50 + j) = (i-1)*100 + j;
        end
    end
    
    % OFF cells: last 50 columns of each row
    ind2 = zeros(50*50, 1);
    for i = 1:50
        for j = 1:50
            ind2((i-1)*50 + j) = (i-1)*100 + (50 + j);
        end
    end
end

function sx = process_fisher_block(block, sx, nspks, fr, F1, F2, ...
                                  noise1, noise2, ind1, ind2, ...
                                  FR, ConPerc, sigma_n, tau_n, dt)
% PROCESS_FISHER_BLOCK - Process one 500ms block for Fisher info
%
% H1 Line: Generate spikes for one stimulus presentation block
%
% Description:
%   Processes a single 500ms block with appropriate stimulus structure
%   for Fisher information analysis.

    Nx = length(fr);
    block_start = (block - 1) * 500;
    
    for t = 1:500
        % Update noise processes
        noise1 = update_noise(noise1, sigma_n, tau_n, dt, length(noise1));
        noise2 = update_noise(noise2, sigma_n, tau_n, dt, length(noise2));
        
        % Calculate firing rates based on stimulus phase
        if t >= 1 && t <= 200
            % Stimulus period - ON cells get baseline + noise
            FR(ind1) = fr(ind1) + F1 * noise1;
            % OFF cells get contrast-modulated rate + noise
            FR(ind2) = fr(ind2) * ConPerc + F2 * noise2;
        elseif t >= 201 && t <= 500
            % Post-stimulus period - very low baseline
            FR(ind1) = 0.005;
            FR(ind2) = 0.005;
        end
        
        % Rectify negative rates
        FR(FR < 0) = 0;
        
        % Generate spikes
        spk = rand(Nx, 1) < FR * dt / 1000;
        
        ns_temp = nnz(spk);
        if ns_temp > 0
            spike_times = block_start + t;
            sx(1, nspks+(1:ns_temp)) = spike_times;
            sx(2, nspks+(1:ns_temp)) = find(spk);
            nspks = nspks + ns_temp;
        end
    end
end

function noise = update_noise(noise, sigma_n, tau_n, dt, N)
% UPDATE_NOISE - Update Ornstein-Uhlenbeck noise process
%
% H1 Line: Update noise vector using Ornstein-Uhlenbeck dynamics
%
% Syntax:
%   noise = update_noise(noise, sigma_n, tau_n, dt, N)
%
% Description:
%   Updates a noise vector following Ornstein-Uhlenbeck dynamics,
%   which produces temporally correlated Gaussian noise.
%
% Inputs:
%   noise   - Current noise state vector
%   sigma_n - Noise standard deviation
%   tau_n   - Noise correlation time constant
%   dt      - Time step
%   N       - Length of noise vector
%
% Outputs:
%   noise - Updated noise vector

    noise = noise + (-noise * dt + sigma_n * randn(N, 1) * sqrt(dt)) / tau_n;
end
