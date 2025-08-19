function [norm_ratio, pdf_density, edges] = get_norm_ind_distr(sim_file, upper_lim, bin_width)
% GET_NORM_IND_DISTR
% Syntax:
%   [norm_ratio, pdf_density, edges] = get_norm_ind_distr(sim_file, upper_lim, bin_width)
% Description:
%   Load three simulations (suffixes _1, _2, _3) and compute the ratio
%   (rate_1 + rate_2) ./ rate_3 per neuron, then return its PDF as a density.
% Inputs:
%   sim_file  : string to a *_3 file; *_1 and *_2 inferred via strrep.
%   upper_lim : positive scalar, max ratio to include in histogram.
%   bin_width : positive scalar, histogram bin width.
% Outputs:
%   norm_ratio : [Ne x 1] vector of ratios; NaNs dropped from PDF calc.
%   pdf_density: row vector, density over [0, upper_lim].
%   edges      : bin edges used.
% Example:
%   [r,p,e] = get_norm_ind_distr('V1fr_Rec_3.mat', 5, 0.1);
% Notes:
%   - Avoids division-by-zero by adding eps to denominator.
%   - Uses param.Ne if present; otherwise infers Ne from data.

  % ---- Validate inputs ----
  arguments
      sim_file (1,1) string
      upper_lim (1,1) double {mustBePositive}
      bin_width (1,1) double {mustBePositive}
  end

  % Infer sibling files
  file1 = strrep(sim_file, 'V1fr_Rec_3', 'V1fr_Rec_1');
  file2 = strrep(sim_file, 'V1fr_Rec_3', 'V1fr_Rec_2');
  file3 = sim_file;

  % Load param and spike trains, compute per-neuron rates (Hz)
  [rate1, Ne, T] = local_load_rate(file1);
  [rate2, ~ ,  ~] = local_load_rate(file2);
  [rate3, ~ ,  ~] = local_load_rate(file3);

  % ---- Ratio ----
  denom = rate3 + eps;                 % avoid divide-by-zero
  norm_ratio = (rate1 + rate2) ./ denom;

  % ---- Histogram of ratio (as density) ----
  valid = norm_ratio < upper_lim & isfinite(norm_ratio);
  edges = 0:bin_width:(upper_lim + eps(upper_lim));
  counts = histcounts(norm_ratio(valid), edges);
  area = sum(counts) * bin_width;
  if area == 0, pdf_density = zeros(size(counts)); else, pdf_density = counts / area; end
end

function [rate, Ne, T] = local_load_rate(f)
  S = load(f, 's2', 'param');
  if ~isfield(S, 's2') || ~isfield(S, 'param')
      error('File %s must contain variables ''s2'' and ''param''.', f);
  end
  if isfield(S.param, 'Ne'), Ne = S.param.Ne; else, Ne = max(S.s2(2,:)); end
  T = S.param.T;
  idx = 1:Ne;
  % Non-overlapping count: Ncount=1, option=1; convert to Hz via /T*1e3 (T in ms)
  counts = spktime2count(S.s2, idx, T, 1, 1);
  rate = counts ./ T * 1e3;
end
