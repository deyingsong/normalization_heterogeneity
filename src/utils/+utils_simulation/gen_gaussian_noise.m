function L = gen_gaussian_noise(T, dt, tau_L, sigma_s)
% GENERATE_GAUSSIAN_NOISE - Generate temporally correlated Gaussian noise
%
% H1 Line: Create Gaussian noise with specified temporal correlation
%
% Syntax:
%   L = generate_gaussian_noise(T, dt, tau_L, sigma_s)
%
% Description:
%   Generates a time series of Gaussian noise with temporal correlations
%   defined by a Gaussian autocovariance function. Uses FFT-based
%   circulant embedding method for efficient generation of long sequences.
%
% Inputs:
%   T       - Total time duration (in time units)
%   dt      - Time step size (in time units)
%   tau_L   - Correlation time constant (in time units)
%   sigma_s - Standard deviation of the noise
%
% Outputs:
%   L - Column vector of noise samples (length = T/dt + 1)
%
% Example:
%   % Generate 10 seconds of noise with 1ms steps and 100ms correlation
%   L = generate_gaussian_noise(10000, 1, 100, 5);
%
% See also: generate_white_noise, generate_ornstein_uhlenbeck_noise

    % Input validation
    validate_noise_parameters(T, dt, tau_L, sigma_s);
    
    % Create time vector
    t = 0:dt:T;
    N = length(t);
    
    % Generate autocovariance function
    C = compute_gaussian_autocovariance(N, dt, tau_L);
    
    % Use circulant embedding for efficiency
    L = generate_noise_circulant_embedding(C, N, sigma_s);
    
    % Ensure column vector output
    L = L(:);
end

function validate_noise_parameters(T, dt, tau_L, sigma_s)
% VALIDATE_NOISE_PARAMETERS - Check noise generation parameters
%
% H1 Line: Validate parameters for Gaussian noise generation
%
% Syntax:
%   validate_noise_parameters(T, dt, tau_L, sigma_s)
%
% Description:
%   Ensures all parameters for noise generation are valid and compatible.
%
% Inputs:
%   T       - Total time (must be positive)
%   dt      - Time step (must be positive and <= T)
%   tau_L   - Correlation time (must be positive)
%   sigma_s - Standard deviation (must be non-negative)

    % Check T
    if ~isscalar(T) || T <= 0
        error('Noise:InvalidT', 'T must be a positive scalar');
    end
    
    % Check dt
    if ~isscalar(dt) || dt <= 0 || dt > T
        error('Noise:InvalidDt', ...
              'dt must be positive and less than T (dt=%.3f, T=%.3f)', dt, T);
    end
    
    % Check tau_L
    if ~isscalar(tau_L) || tau_L <= 0
        error('Noise:InvalidTau', 'tau_L must be a positive scalar');
    end
    
    % Check sigma_s
    if ~isscalar(sigma_s) || sigma_s < 0
        error('Noise:InvalidSigma', 'sigma_s must be non-negative');
    end
    
    % Warning for potential issues
    if tau_L < dt
        warning('Noise:ShortCorrelation', ...
                'Correlation time (%.3f) is less than time step (%.3f)', ...
                tau_L, dt);
    end
    
    if T/dt > 1e6
        warning('Noise:LargeArray', ...
                'Generating %.1e samples - may require significant memory', T/dt);
    end
end

function C = compute_gaussian_autocovariance(N, dt, tau_L)
% COMPUTE_GAUSSIAN_AUTOCOVARIANCE - Calculate Gaussian autocovariance
%
% H1 Line: Compute autocovariance function for Gaussian process
%
% Syntax:
%   C = compute_gaussian_autocovariance(N, dt, tau_L)
%
% Description:
%   Computes the autocovariance function C(tau) = exp(-tau^2/(2*tau_L^2))
%   for generating temporally correlated noise.
%
% Inputs:
%   N     - Number of time points
%   dt    - Time step
%   tau_L - Correlation time constant
%
% Outputs:
%   C - Autocovariance vector

    tau = (-N+1:N-1) * dt;
    C = exp(-tau.^2 / (2 * tau_L^2));
end

function L = generate_noise_circulant_embedding(C, N, sigma_s)
% GENERATE_NOISE_CIRCULANT_EMBEDDING - Generate noise using FFT method
%
% H1 Line: Use circulant embedding to generate correlated noise
%
% Syntax:
%   L = generate_noise_circulant_embedding(C, N, sigma_s)
%
% Description:
%   Efficiently generates correlated Gaussian noise using the
%   circulant embedding method with FFT.
%
% Inputs:
%   C       - Autocovariance vector
%   N       - Number of samples needed
%   sigma_s - Standard deviation
%
% Outputs:
%   L - Noise time series

    % Create circulant covariance matrix representation
    C_fft = [C(N:end), 0, C(1:N-1)];  % Make symmetric around 0
    
    % Compute spectral density via FFT
    S = real(fft(C_fft));
    
    % Check for negative eigenvalues (shouldn't happen for valid covariance)
    if any(S < -1e-10)
        error('Noise:NegativeEigenvalues', ...
              'Covariance matrix is not positive semi-definite');
    end
    
    % Small negative values due to numerical errors - set to zero
    S(S < 0) = 0;
    
    % Generate white Gaussian noise
    W = randn(1, length(S));
    
    % Filter in frequency domain
    L_freq = sqrt(S) .* fft(W);
    
    % Transform back to time domain
    L = real(ifft(L_freq));
    
    % Extract relevant samples and scale
    L = L(1:N) * sigma_s;
end
