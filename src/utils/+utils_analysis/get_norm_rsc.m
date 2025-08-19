function [dataheat, ind_lim] = get_norm_rsc(spktr, norm_ind, nbins, prct_bounds)
% GET_NORM_RSC
% Syntax:
%   [dataheat, ind_lim] = get_norm_rsc(spktr, norm_ind)
%   [dataheat, ind_lim] = get_norm_rsc(spktr, norm_ind, nbins, prct_bounds)
% Description:
%   Given neurons-by-time spike counts (spktr) and a per-neuron normalization
%   index (norm_ind), compute a heat map of mean pairwise correlations across
%   bins of the index.
% Inputs:
%   spktr       : [N x T] or [T x N] matrix of spike counts per bin.
%   norm_ind    : [N x 1] vector of index values for each neuron.
%   nbins       : (optional) number of bins (default 20).
%   prct_bounds : (optional) [low high] percentiles to define range (default [5 95]).
% Outputs:
%   dataheat : [nbins x nbins] mean correlation per (i,j) bin pair.
%   ind_lim  : [1x2] bounds used for binning.
% Example:
%   H = get_norm_rsc(spktr, ratio);

  if nargin < 3 || isempty(nbins), nbins = 20; end
  if nargin < 4 || isempty(prct_bounds), prct_bounds = [5 95]; end
  validateattributes(nbins, {'numeric'}, {'scalar','integer','>=',2});
  validateattributes(prct_bounds, {'numeric'}, {'vector','numel',2,'increasing'});

  % Ensure [N x T]
  if size(spktr,1) ~= numel(norm_ind)
      spktr = spktr';
      if size(spktr,1) ~= numel(norm_ind)
          error('%s: spktr shape does not match norm_ind length.', mfilename);
      end
  end

  low  = prctile(norm_ind, prct_bounds(1), 'all');
  high = prctile(norm_ind, prct_bounds(2), 'all');
  ind_lim = [low, high];
  edges = linspace(low, high, nbins+1);

  dataheat = zeros(nbins, nbins);

  for ii = 1:nbins
      mask_i = norm_ind >= edges(ii) & norm_ind < edges(ii+1);
      Xi = spktr(mask_i, :);
      for jj = 1:nbins
          mask_j = norm_ind >= edges(jj) & norm_ind < edges(jj+1);
          Xj = spktr(mask_j, :);
          if isempty(Xi) || isempty(Xj)
              dataheat(ii,jj) = NaN; continue; %#ok<AGROW>
          end
          if ii == jj
              C = cov(Xi'); R = corrcov(C); r = R(triu(true(size(R)),1));
              dataheat(ii,jj) = mean(r, 'omitnan');
          else
              R = corr(Xi', Xj');
              dataheat(ii,jj) = mean(R(:), 'omitnan');
          end
      end
  end
end
