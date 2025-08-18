function [fr_pdf, fr_mean, fr_std, edges] = get_fr_pdf_mean_std(fr, fr_max, bin_width)
% GET_FR_PDF_MEAN_STD
% Syntax:
%   [fr_pdf, fr_mean, fr_std, edges] = get_fr_pdf_mean_std(fr, fr_max, bin_width)
% Description:
%   Compute the firing-rate PDF (normalized as a density), mean, and std.
%   Avoids magic numbers by making the max value and bin width explicit.
% Inputs:
%   fr        : numeric vector of firing rates (Hz). NaNs ignored.
%   fr_max    : positive scalar, right edge (inclusive) for histogram.
%   bin_width : positive scalar, width of histogram bins (Hz).
% Outputs:
%   fr_pdf : row vector, probability density per bin (area sums to 1).
%   fr_mean: scalar, mean firing rate (ignoring NaNs).
%   fr_std : scalar, std of firing rate (ignoring NaNs).
%   edges  : vector, bin edges used by histcounts.
% Example:
%   fr = [2.1 5 3 7 8 NaN];
%   [pdf, m, s, edges] = get_fr_pdf_mean_std(fr, 10, 1);
% Notes (performance):
%   - For profiling, wrap calls with tic/toc or use PROFILE ON/OFF.
%   - Avoids copying large arrays; operates directly on inputs.

  % ---- Input validation ----
  validateattributes(fr, {'numeric'}, {'vector'}, mfilename, 'fr', 1);
  validateattributes(fr_max, {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'fr_max', 2);
  validateattributes(bin_width, {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'bin_width', 3);

  fr = fr(:)';                         % ensure row vector
  fr = fr(~isnan(fr));                  % drop NaNs (documented)
  if isempty(fr)
      fr_pdf = []; fr_mean = NaN; fr_std = NaN; edges = 0:bin_width:fr_max; return;
  end

  % ---- Histogram as density ----
  % Use inclusive rightmost edge by adding eps.
  edges = 0:bin_width:(fr_max + eps(fr_max));
  counts = histcounts(fr, edges);
  area = sum(counts) * bin_width;
  if area == 0
      fr_pdf = zeros(size(counts));
  else
      fr_pdf = counts / area;          % density so integral = 1
  end

  % ---- Summary statistics ----
  fr_mean = mean(fr, 'omitnan');
  fr_std  = std(fr, 0, 'omitnan');
end
