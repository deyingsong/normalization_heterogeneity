function x = gen_gamma_integer_array(alpha, theta, L, N)
% GENERATE_GAMMA_INTEGER_ARRAY - Generate integers from gamma distribution
%
% H1 Line: Generate L integers summing to N from gamma distribution
%
% Syntax:
%   x = generate_gamma_integer_array(alpha, theta, L, N)
%
% Description:
%   Generates an array of L non-negative integers that sum to exactly N,
%   where the relative proportions follow a gamma distribution with shape
%   parameter alpha and scale parameter theta. Uses intelligent rounding
%   to preserve the total sum while maintaining distribution properties.
%
% Inputs:
%   alpha - Shape parameter of gamma distribution (alpha > 0)
%   theta - Scale parameter of gamma distribution (theta > 0)
%   L     - Number of integers to generate
%   N     - Target sum of all integers
%
% Outputs:
%   x - Row vector of L non-negative integers summing to N
%
% Example:
%   % Generate 10 integers summing to 100 from gamma(2, 3)
%   x = generate_gamma_integer_array(2, 3, 10, 100);
%   
% See also: gamrnd, round_preserving_sum

    % Input validation
    validate_gamma_parameters(alpha, theta, L, N);
    
    % Generate samples from gamma distribution
    gamma_samples = gamrnd(alpha, theta, [1, L]);
    
    % Handle edge case of all zeros (extremely rare but possible)
    if sum(gamma_samples) == 0
        warning('GammaInt:AllZeros', ...
                'All gamma samples were zero, using uniform distribution');
        gamma_samples = ones(1, L);
    end
    
    % Normalize to create probability distribution
    p = gamma_samples / sum(gamma_samples);
    
    % Round preserving the sum
    x = round_preserving_sum(p * N, N);
end

function validate_gamma_parameters(alpha, theta, L, N)
% VALIDATE_GAMMA_PARAMETERS - Validate inputs for gamma integer generation
%
% H1 Line: Check validity of gamma distribution parameters
%
% Syntax:
%   validate_gamma_parameters(alpha, theta, L, N)
%
% Description:
%   Validates that all parameters for gamma integer array generation
%   are within acceptable ranges and have compatible types.
%
% Inputs:
%   alpha - Shape parameter (must be positive)
%   theta - Scale parameter (must be positive)
%   L     - Array length (must be positive integer)
%   N     - Target sum (must be non-negative integer)
%
% Example:
%   validate_gamma_parameters(2, 3, 10, 100);

    % Check alpha
    if ~isscalar(alpha) || ~isnumeric(alpha) || alpha <= 0
        error('GammaInt:InvalidAlpha', ...
              'Alpha must be a positive scalar (got alpha=%.2f)', alpha);
    end
    
    % Check theta
    if ~isscalar(theta) || ~isnumeric(theta) || theta <= 0
        error('GammaInt:InvalidTheta', ...
              'Theta must be a positive scalar (got theta=%.2f)', theta);
    end
    
    % Check L
    if ~isscalar(L) || ~isnumeric(L) || L <= 0 || floor(L) ~= L
        error('GammaInt:InvalidL', ...
              'L must be a positive integer (got L=%.2f)', L);
    end
    
    % Check N
    if ~isscalar(N) || ~isnumeric(N) || N < 0 || floor(N) ~= N
        error('GammaInt:InvalidN', ...
              'N must be a non-negative integer (got N=%.2f)', N);
    end
    
    % Check compatibility
    if N > 0 && L > N
        warning('GammaInt:SparseSolution', ...
                'L > N will result in many zeros (L=%d, N=%d)', L, N);
    end
end

function x = round_preserving_sum(x_real, target_sum)
% ROUND_PRESERVING_SUM - Round array while preserving exact sum
%
% H1 Line: Round real numbers to integers preserving total sum
%
% Syntax:
%   x = round_preserving_sum(x_real, target_sum)
%
% Description:
%   Rounds an array of real numbers to integers such that the sum
%   equals exactly target_sum. Uses largest remainder method with
%   random tie-breaking for fairness.
%
% Inputs:
%   x_real     - Array of real numbers to round
%   target_sum - Desired sum of rounded integers
%
% Outputs:
%   x - Array of integers with same size as x_real
%
% Example:
%   x = round_preserving_sum([2.3, 4.7, 1.8], 9);

    % Initial rounding (floor)
    x = floor(x_real);
    remainder = target_sum - sum(x);
    
    % Handle negative remainder (should not happen with valid inputs)
    if remainder < 0
        error('RoundSum:NegativeRemainder', ...
              'Cannot achieve target sum - check inputs');
    end
    
    % Distribute remainder to elements with largest fractional parts
    if remainder > 0
        fractional_parts = x_real - x;
        
        % Add small random values to break ties fairly
        tie_breaker = rand(size(fractional_parts)) * 1e-10;
        fractional_parts = fractional_parts + tie_breaker;
        
        % Sort and increment top 'remainder' elements
        [~, idx] = sort(fractional_parts, 'descend');
        
        % Ensure we don't try to increment more elements than exist
        remainder = min(remainder, length(x));
        x(idx(1:remainder)) = x(idx(1:remainder)) + 1;
    end
    
    % Verify sum (for debugging)
    if sum(x) ~= target_sum
        error('RoundSum:SumMismatch', ...
              'Sum mismatch after rounding: got %d, expected %d', ...
              sum(x), target_sum);
    end
end
