function [dataheat, ind_lim] = get_norm_rsc_simfile(sim_file, ind_lim)
% GET_NORM_RSC_SIMFILE
% Syntax:
%   [dataheat, ind_lim] = get_norm_rsc_simfile(sim_file, ind_lim)
% Description:
%   From a simulation file and its *_1/*_2 siblings, compute per-neuron
%   ratio (rate_1+rate_2)/rate_3, select neurons within percentile bounds,
%   smooth their time-binned spike counts, and compute mean pairwise
%   correlations across ratio bins to produce a heat map.
% Inputs:
%   sim_file : string to a *_3 file; *_1 and *_2 inferred via strrep.
%   ind_lim  : [] or [low high] ratio bounds. If empty, uses [5,95] percentiles.
% Outputs:
%   dataheat : [nbins x nbins] matrix of mean correlations.
%   ind_lim  : [1x2] bounds actually used.
% Example:
%   [H, lim] = get_norm_rsc_simfile('V1fr_Rec_3.mat', []);
% Parameters (constants made explicit):
%   NBINS             = 20;            % number of ratio bins
%   BIN_MS            = 50;            % spike-count bin size (ms)
%   SMOOTH_TAPS      = 4;              % moving-average taps
%   PERCENTILE_RANGE = [5, 95];        % default ratio bounds

  % ---- Constants ----
  NBINS = 20; BIN_MS = 50; SMOOTH_TAPS = 4; PERCENTILE_RANGE = [5, 95];

  % ---- Load spike trains and build smoothed spike count matrix ----
  S3 = load(sim_file, 's2', 'param');
  if ~isfield(S3, 's2') || ~isfield(S3, 'param'), error('s2/param missing in %s', sim_file); end
  Ne = getfieldwithdefault(S3.param, 'Ne', max(S3.s2(2,:)));
  T  = S3.param.T;

  Nt = floor(T / BIN_MS);                     % #bins (non-overlapping)
  counts = spktime2count(S3.s2, 1:Ne, BIN_MS, Nt, 1);   % [Ne x Nt]

  % Discard the first 1000/BIN_MS bins if T is in ms and first second is transient
  transient_bins = max(0, floor(1000 / BIN_MS));
  if Nt > transient_bins
      counts = counts(:, transient_bins+1:end);
  end

  % Smooth across time with a simple moving average of length SMOOTH_TAPS
  if SMOOTH_TAPS > 1
      h = ones(1, SMOOTH_TAPS) / SMOOTH_TAPS;
      counts = conv2(counts, h, 'same');
  end

  % ---- Compute per-neuron rates (Hz) in *_1, *_2, *_3 ----
  [rate1, ~] = local_rate_from_file(strrep(sim_file, 'V1fr_Rec_3', 'V1fr_Rec_1'));
  [rate2, ~] = local_rate_from_file(strrep(sim_file, 'V1fr_Rec_3', 'V1fr_Rec_2'));
  [rate3, ~] = local_rate_from_file(sim_file);

  ratio = (rate1 + rate2) ./ (rate3 + eps);

  % ---- Subset neurons within ratio bounds ----
  if isempty(ind_lim)
      low  = prctile(ratio, PERCENTILE_RANGE(1), 'all');
      high = prctile(ratio, PERCENTILE_RANGE(2), 'all');
      ind_lim = [low, high];
  else
      validateattributes(ind_lim, {'numeric'}, {'vector','numel',2,'increasing'});
      low = ind_lim(1); high = ind_lim(2);
  end

  mask = ratio >= low & ratio < high & all(isfinite([rate1, rate2, rate3]), 2);
  spktr = counts(mask, :);
  norm_vals = ratio(mask);

  % Ensure rows correspond to neurons
  if size(spktr,1) ~= numel(norm_vals)
      spktr = spktr';
      if size(spktr,1) ~= numel(norm_vals)
          error('%s: spike/time matrix shape mismatch after transpose.', mfilename);
      end
  end

  % ---- Heat map of mean correlations across ratio bins ----
  edges = linspace(low, high, NBINS+1);
  dataheat = zeros(NBINS, NBINS);

  for ii = 1:NBINS
      sel_i = norm_vals >= edges(ii) & norm_vals < edges(ii+1);
      Xi = spktr(sel_i, :);
      for jj = 1:NBINS
          sel_j = norm_vals >= edges(jj) & norm_vals < edges(jj+1);
          Xj = spktr(sel_j, :);
          if isempty(Xi) || isempty(Xj)
              dataheat(ii, jj) = NaN;  %#ok<AGROW>
          elseif ii == jj
              C = cov(Xi');  R = corrcov(C);  rr = R(triu(true(size(R)),1));
              dataheat(ii, jj) = mean(rr, 'omitnan');
          else
              R = corr(Xi', Xj');
              dataheat(ii, jj) = mean(R(:), 'omitnan');
          end
      end
  end
end

function [rate, T] = local_rate_from_file(f)
  S = load(f, 's2', 'param');
  if ~isfield(S, 's2') || ~isfield(S, 'param')
      error('File %s must contain variables ''s2'' and ''param''.', f);
  end
  Ne = getfieldwithdefault(S.param, 'Ne', max(S.s2(2,:)));
  T  = S.param.T;
  rate = spktime2count(S.s2, 1:Ne, T, 1, 1) ./ T * 1e3; % Hz
end

function v = getfieldwithdefault(S, name, default)
  if isfield(S, name), v = S.(name); else, v = default; end
end
