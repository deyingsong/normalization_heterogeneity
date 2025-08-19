function theta = generate_ori_map(Nx, Kc, varargin)
% GENERATE_ORI_MAP  Orientation map generator (Kaschube et al., 2010, Eq. 20).
%
% Syntax
%   theta = generate_ori_map(Nx, Kc, Name,Value, ...)
%
% Description
%   Implements the orientation map synthesis from Kaschube et al. (2010) Supp.
%   Materials, Eq. 20. Produces a square map of size [Nx x Nx] with preferred
%   orientations in [0,1] (i.e., fraction of pi).
%
% Inputs
%   Nx : Positive integer. Map side length (Nx x Nx grid).
%   Kc : Positive real. Spatial frequency parameter (|k| = Kc).
%
% Name-Value Pairs (optional)
%   'NumWaves'  : Positive integer, number of plane waves (default 30).
%   'Seed'      : Numeric scalar. If provided, rng(Seed,'twister') is set for
%                 reproducibility (affects random signs and phases).
%   'ShowPlot'  : Logical. If true, displays the map (default false).
%
% Output
%   theta : [Nx x Nx] double, preferred orientation \in [0,1].
%
% Example
%   theta = generate_ori_map(50, 5*2*pi, 'NumWaves', 30, 'Seed', 1, 'ShowPlot', true);
%
% Author: Refactored for clarity & robustness.
% -------------------------------------------------------------------------

% Parse inputs
p = inputParser;
p.FunctionName = mfilename;
addRequired(p,'Nx', @(x)isnumeric(x)&&isscalar(x)&&x>0&&mod(x,1)==0);
addRequired(p,'Kc', @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'NumWaves', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0&&mod(x,1)==0);
addParameter(p,'Seed', [], @(x)isnumeric(x)&&isscalar(x) || isempty(x));
addParameter(p,'ShowPlot', false, @(x)islogical(x)&&isscalar(x));
parse(p,Nx,Kc,varargin{:});
S = p.Results;

if ~isempty(S.Seed)
    rng(S.Seed,'twister');
end

Ny = Nx;
dx = 1/Nx;
x = 0:dx:1;  % length Nx+1; match original intent using lattice including both ends
y = x;

n = S.NumWaves;
% k-vectors uniformly spaced directions with magnitude Kc
angles = (0:n-1).' * pi/n;
k = Kc * [cos(angles), sin(angles)];

% Random signs l \in {-1,+1} (no zeros)
l = sign(rand(n,1)-0.5);
l(l==0) = 1;

% Random phases phi \in [0, 2*pi)
phi = rand(n,1) * 2*pi;

% Build complex sum z at each lattice point
z = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        pos = [x(i); y(j)];
        z(i,j) = sum( exp(1i * ( l .* (k*pos) + phi )) );
    end
end

% Preferred orientation in [0,1] (fraction of pi)
theta = angle(z) / pi / 2 + 0.5;

if any(isnan(theta(:)))
    error('%s:NaNOutput','Output contains NaN values. Please check inputs.');
end

if S.ShowPlot
    figure('Name','Orientation map','Color','w');
    imagesc(theta);
    axis image off;
    colormap(hsv);
    colorbar;
    title('Preferred Orientation (fraction of \pi)');
end
