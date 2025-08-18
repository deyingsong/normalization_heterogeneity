function [rsc_mean, rsc_std] = get_rsc_mean_std(spktr)
% GET_RSC_MEAN_STD
% Syntax:
%   [rsc_mean, rsc_std] = get_rsc_mean_std(spktr)
% Description:
%   Compute mean and standard deviation of pairwise correlations among
%   neurons given spike-count matrix.
% Input:
%   spktr : [N x T] spike-count matrix (rows = neurons). NaNs tolerated.
% Outputs:
%   rsc_mean : scalar mean of off-diagonal correlations.
%   rsc_std  : scalar std  of off-diagonal correlations.
% Example:
%   [m,s] = get_rsc_mean_std(spktr);

  C = cov(spktr');
  R = corrcov(C);
  r = R(triu(true(size(R)),1));
  r = r(~isnan(r));
  rsc_mean = mean(r);
  rsc_std  = std(r);
end
