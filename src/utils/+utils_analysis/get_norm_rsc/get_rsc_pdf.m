function [rsc_pdf, edges] = get_rsc_pdf(spktr, bin_width)
% GET_RSC_PDF
% Syntax:
%   [rsc_pdf, edges] = get_rsc_pdf(spktr)
%   [rsc_pdf, edges] = get_rsc_pdf(spktr, bin_width)
% Description:
%   Return the probability density of pairwise correlations across neurons.
% Inputs:
%   spktr     : [N x T] spike-count matrix.
%   bin_width : optional scalar (default 0.04) for histogram bin width.
% Outputs:
%   rsc_pdf : row vector density across [-1, 1].
%   edges   : histogram edges used.
% Example:
%   [pdf, e] = get_rsc_pdf(spktr, 0.05);

  if nargin < 2 || isempty(bin_width), bin_width = 0.04; end
  validateattributes(bin_width, {'numeric'}, {'scalar','positive','<=',1});

  C = cov(spktr');
  R = corrcov(C);
  r = R(triu(true(size(R)),1));

  edges = -1:bin_width:(1 + eps(1));
  counts = histcounts(r, edges);
  area = sum(counts) * bin_width;
  if area == 0
      rsc_pdf = zeros(size(counts));
  else
      rsc_pdf = counts / area;
  end
end

