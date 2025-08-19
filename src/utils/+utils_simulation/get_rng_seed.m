function rng_seed = get_rng_seed(options)
% GET_RNG_SEED - Get or generate random seed
%
% H1 Line: Obtain random seed for reproducibility
%
% Syntax:
%   rng_seed = get_rng_seed(options)
%
% Description:
%   Returns the random seed from options or generates a default one.
%
% Inputs:
%   options - Options structure
%
% Outputs:
%   rng_seed - Random number generator seed
%
% Example:
%   seed = get_rng_seed(options);

    if isfield(options, 'rng_seed')
        rng_seed = options.rng_seed;
    else
        rng_seed = 80001; % Default seed for reproducibility
    end
end
