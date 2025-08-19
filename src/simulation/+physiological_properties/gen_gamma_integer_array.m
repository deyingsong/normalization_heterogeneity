function x = gen_gamma_integer_array(alpha, theta, L, N)
% ATTTREC_GENERATE_GAMMA_INTEGER_ARRAY  Allocate N integers via Gamma-proportional rounding.
%
% Syntax
%   x = atttRec_generate_gamma_integer_array(alpha, theta, L, N)
%
% Description
%   Draw L i.i.d. Gamma(shape=alpha, scale=theta) samples, convert them into a
%   probability vector p, and allocate a total integer count N by flooring N*p
%   and distributing any remainder to the largest fractional parts (ties broken
%   randomly). This preserves sum(x) = N and produces a 1xL nonnegative integer row.
%
% Inputs
%   alpha : Positive scalar (Gamma shape).
%   theta : Positive scalar (Gamma scale).
%   L     : Positive integer. Number of categories/samples.
%   N     : Nonnegative integer. Total count to allocate.
%
% Output
%   x     : 1xL row vector of nonnegative integers with sum(x) == N.
%
% Example
%   rng(1,'twister');                           % for reproducibility (optional)
%   x = atttRec_generate_gamma_integer_array(2, 0.5, 10, 123)
%
% Notes
%   * Requires Statistics and Machine Learning Toolbox for GAMRND.
%   * If GAMRND is unavailable or the sampled mass is numerically zero,
%     the function falls back to a uniform allocation baseline.
%
% See also  gamrnd, rng
% -------------------------------------------------------------------------

% ---------- validation ----------
validateattributes(alpha, {'double','single'},{'real','finite','scalar','>',0}, mfilename, 'alpha');
validateattributes(theta, {'double','single'},{'real','finite','scalar','>',0}, mfilename, 'theta');
validateattributes(L,     {'double','single'},{'real','finite','scalar','integer','>=',1}, mfilename, 'L');
validateattributes(N,     {'double','single'},{'real','finite','scalar','integer','>=',0}, mfilename, 'N');

% ---------- parameters ----------
TIE_EPS = 1e-9;  % small jitter for random tie-breaking (no "magic": named & documented)

% ---------- sample from Gamma(alpha, theta) ----------
gamma_samples = [];
try
    gamma_samples = gamrnd(alpha, theta, [1, L]);
catch ME
    if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
        error('%s:MissingGamrnd', ['gamrnd is unavailable. This function requires the ' ...
            'Statistics and Machine Learning Toolbox. Install it, or replace the sampler.']);
    else
        rethrow(ME);
    end
end

% Guard against numerical pathologies
if ~isfloat(gamma_samples) || any(~isfinite(gamma_samples))
    error('%s:InvalidSamples','Gamma sampler returned NaN/Inf values. Check alpha/theta.', mfilename);
end
if all(gamma_samples <= 0)
    % fallback to uniform mass if degenerately zero (extremely unlikely for valid alpha/theta)
    gamma_samples = ones(1, L);
end

% ---------- convert to probabilities ----------
total_mass = sum(gamma_samples);
if ~isfinite(total_mass) || total_mass <= 0
    % uniform fallback
    p = ones(1, L) / L;
else
    p = gamma_samples ./ total_mass;
end

% ---------- integer allocation ----------
x_real   = N * p;                 % target real-valued counts
x        = floor(x_real);         % base allocation
remainder = N - sum(x);           % how many integers still to place (should be in [0, L])

% Numerical safety: clamp remainder to [0, L]
remainder = max(0, min(L, round(double(remainder))));

if remainder > 0
    frac = x_real - x;                            % fractional parts in [0,1)
    % Randomized tie-breaking so equal fractions are split stochastically
    [~, idx] = sort(frac + rand(1, L) * TIE_EPS, 'descend');
    x(idx(1:remainder)) = x(idx(1:remainder)) + 1;
end

% ---------- final checks & shape ----------
if sum(x) ~= N
    % Very defensive (should not happen). Repair with minimal adjustments.
    delta = N - sum(x);
    if delta > 0
        % add ones to largest p entries
        [~, idx] = sort(p, 'descend');
        k = 1;
        while delta > 0
            x(idx(k)) = x(idx(k)) + 1;
            delta = delta - 1;
            k = k + 1; if k > L, k = 1; end
        end
    elseif delta < 0
        % remove ones from largest entries
        [~, idx] = sort(x, 'descend');
        k = 1;
        while delta < 0 && any(x)
            if x(idx(k)) > 0
                x(idx(k)) = x(idx(k)) - 1;
                delta = delta + 1;
            end
            k = k + 1; if k > L, k = 1; end
        end
    end
end

x = reshape(x, 1, L);  % ensure row vector

end
