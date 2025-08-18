function rsc_dist = get_rsc_dist(spktr, neuron_indices, dist_max, grid_size, bin_subdivisions)
% GET_RSC_DIST
% Syntax:
%   rsc_dist = get_rsc_dist(spktr, neuron_indices, dist_max)
%   rsc_dist = get_rsc_dist(spktr, neuron_indices, dist_max, grid_size, bin_subdivisions)
% Description:
%   Compute mean pairwise spike-count correlation (rSC) as a function of
%   toroidal distance between neuron positions on a grid.
% Inputs:
%   spktr            : [N x T] spike-count matrix (rows = neurons).
%   neuron_indices   : [N x 1] linear indices of neuron positions on a grid.
%   dist_max         : positive scalar, max grid distance to evaluate.
%   grid_size        : (optional) [W H] grid dimensions (default [200 100]).
%   bin_subdivisions : (optional) integer, #sub-bins per distance unit (default 100).
% Output:
%   rsc_dist : [dist_max*bin_subdivisions x 1] mean rSC per fine distance bin.
% Example:
%   r = get_rsc_dist(spktr, ind1, 80, [200 100], 100);

  if nargin < 4 || isempty(grid_size),        grid_size = [200 100]; end
  if nargin < 5 || isempty(bin_subdivisions), bin_subdivisions = 100; end
  validateattributes(dist_max, {'numeric'}, {'scalar','positive'});
  validateattributes(grid_size, {'numeric'}, {'vector','numel',2,'integer','>=',1});
  validateattributes(bin_subdivisions, {'numeric'}, {'scalar','integer','>=',1});

  % Pairwise correlation from spike counts
  C = cov(spktr');
  R = corrcov(C);
  R = R(triu(true(size(R)),1));        % vector of upper-triangle correlations

  % Toroidal distances between neuron indices
  W = grid_size(1); H = grid_size(2);
  x = mod(neuron_indices-1, W) + 1;    % 1..W
  y = ceil(neuron_indices / W);        % 1..H
  dx = min(abs(x - x'), W - abs(x - x'));
  dy = min(abs(y - y'), H - abs(y - y'));
  D  = sqrt(dx.^2 + dy.^2);
  D  = D(triu(true(size(D)),1));       % vectorize upper triangle

  % Fine distance histogram and mean rSC per bin
  nb = dist_max * bin_subdivisions;
  rsc_dist = zeros(nb, 1);
  for i = 1:nb
      lo = (i-1) / bin_subdivisions;
      hi = i / bin_subdivisions;
      sel = (D > lo) & (D <= hi);
      rsc_dist(i) = mean(R(sel), 'omitnan');
  end
end
