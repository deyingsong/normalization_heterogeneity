function [F, F1, F2] = generate_Gabor_filter(ori_map_V11, ori_map_V12, varargin)
% GENERATE_GABOR_FILTER  Build population Gabor filters from two orientation maps.
% 
% Syntax
%   [F, F1, F2] = generate_Gabor_filter(ori_map_V11, ori_map_V12, Name,Value, ...)
%
% Description
%   This function constructs bank(s) of Gabor-like filters (F1,F2) whose
%   preferred orientations are taken from two orientation maps (e.g., two V1
%   sheets). It then stacks/tiles them into a single matrix F matching the
%   original author's layout (two sheets tiled along rows, two sub-banks
%   tiled along columns).
%
%   Both ori_map_V11 and ori_map_V12 can be either:
%     * a numeric matrix of size [Ny, Nx] containing preferred orientations
%       in [0,1] (fraction of pi), or
%     * a string/char path to a .mat file. By default, the function will
%       look for variables 'theta_map1' and 'theta_map2' respectively; you
%       may override the variable names via Name-Value pairs.
%
% Inputs
%   ori_map_V11 : numeric matrix or char/string (path to .mat). Orientation map #1.
%   ori_map_V12 : numeric matrix or char/string (path to .mat). Orientation map #2.
%
% Name-Value Pairs (all optional)
%   'Sigma'     : Positive scalar. Gaussian envelope std of the Gabor (default 0.2).
%   'Lambda'    : Positive scalar. Wavelength of the cosine carrier (default 0.6).
%   'Dx'        : Spatial sampling step for the Gabor patch (default 0.07).
%   'Extent'    : Half-extent of the square patch in each axis (default 0.49).
%                 The patch spans [-Extent : Dx : +Extent].
%   'Theta1'    : Scalar (radians). Test angle used only for internal
%                 normalization of F1 (default pi).
%   'Theta2'    : Scalar (radians). Test angle used only for internal
%                 normalization of F2 (default pi/2).
%   'rX'        : Positive scalar scaling factor used by the original code
%                 to normalize average response; retained for compatibility
%                 (default 0.01).
%   'Seed'      : Numeric scalar. If provided, rng(Seed,'twister') is set
%                 to ensure reproducibility of any randomized operations.
%                 (This function is otherwise deterministic.)
%   'VarName1'  : Name of the variable to load from ori_map_V11 if it is a
%                 .mat file (default 'theta_map1').
%   'VarName2'  : Name of the variable to load from ori_map_V12 if it is a
%                 .mat file (default 'theta_map2').
%
% Outputs
%   F  : [2*Ny*Ny  , 2*K] double. Tiled bank of filters following the original
%        layout (two cortical sheets stacked by rows, two sub-banks by columns).
%   F1 : [Ny*Ny    , K  ] double. Filter bank derived from ori_map_V11.
%   F2 : [Ny*Ny    , K  ] double. Filter bank derived from ori_map_V12.
%
% Example
%   % Example using .mat files with variables 'theta_map1' and 'theta_map2'
%   % and default parameters. (Set the random seed for reproducibility.)
%   [F, F1, F2] = generate_Gabor_filter('V11.mat','V12.mat', ...
%                     'Seed', 1, 'Sigma', 0.2, 'Lambda', 0.6, ...
%                     'Dx', 0.07, 'Extent', 0.49);
%
% Author: Refactored for clarity & robustness.
% -------------------------------------------------------------------------

% ---- Parse inputs
p = inputParser;
p.FunctionName = mfilename;
addRequired(p,'ori_map_V11',@(x)isnumeric(x) || ischar(x) || isstring(x));
addRequired(p,'ori_map_V12',@(x)isnumeric(x) || ischar(x) || isstring(x));
addParameter(p,'Sigma',  0.2,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Lambda', 0.6,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Dx',     0.07,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Extent', 0.49,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Theta1', pi,    @(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Theta2', 0.5*pi,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'rX',     0.01,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Seed',   [],    @(x)isnumeric(x)&&isscalar(x) || isempty(x));
addParameter(p,'VarName1','theta_map1',@(s)ischar(s)||isstring(s));
addParameter(p,'VarName2','theta_map2',@(s)ischar(s)||isstring(s));
parse(p,ori_map_V11,ori_map_V12,varargin{:});
S = p.Results;

% ---- Reproducibility (no randomness used here, but keep for consistency)
if ~isempty(S.Seed)
    rng(S.Seed,'twister');
end

% ---- Resolve/validate orientation maps
theta_map1 = resolveOrientationMap(ori_map_V11, S.VarName1);
theta_map2 = resolveOrientationMap(ori_map_V12, S.VarName2);

validateattributes(theta_map1, {'double','single'}, {'2d','nonempty','>=',0,'<=',1}, mfilename, 'theta_map1');
validateattributes(theta_map2, {'double','single'}, {'2d','nonempty','>=',0,'<=',1}, mfilename, 'theta_map2');

if any(isnan(theta_map1(:))) || any(isnan(theta_map2(:)))
    error('%s:NaNInput','Orientation maps contain NaN values. Please clean or interpolate them before use.');
end

[Ny1, Nx1] = size(theta_map1);
[Ny2, Nx2] = size(theta_map2);
if Ny1~=Nx1 || Ny2~=Nx2 || Ny1~=Ny2 || Nx1~=Nx2
    error('%s:MapDims','Expected both orientation maps to be square and of equal size. Got [%dx%d] and [%dx%d].', ...
        mfilename, Ny1, Nx1, Ny2, Nx2);
end
Ny = Ny1; % map side length

% ---- Build Gabor patch grid
x = -S.Extent:S.Dx:S.Extent;  % e.g., -0.49:0.07:0.49  => length Kside
[X, Y] = meshgrid(x, x);
Xv = X(:); Yv = Y(:);
K = numel(Xv); % number of pixels per patch

% ---- Define Gabor responses as anonymous functions (vectorized over theta)
Imag   = @(theta) exp(-(Xv.^2+Yv.^2)/(2*S.Sigma^2)) .* cos( 2*pi/S.Lambda*(Xv.*cos(theta) + Yv.*sin(theta)) );
% Normalize each filter by energy of a reference (theta=0) to keep scale comparable
normDen = sum(Imag(0).^2);
Filter  = @(theta) (exp(-(Xv.^2+Yv.^2)/(2*S.Sigma^2)) .* cos( 2*pi/S.Lambda*(Xv.*cos(theta) + Yv.*sin(theta)) )) ./ normDen;

% ---- Flatten orientation maps into [Ny^2 x 1] angle (radians)
theta_vec1 = reshape(theta_map1, [], 1) * pi;
theta_vec2 = reshape(theta_map2, [], 1) * pi;

% ---- Build filter banks (rows = neurons/pixels of map; cols = patch pixels)
F1 = Filter(theta_vec1.')';  % [Ny^2 x K]
fr = F1 * Imag(S.Theta1);    % mean response at test angle
scale = S.rX ./ mean(fr);    % single scalar scale used in original code
F1 = F1 * scale;

F2 = Filter(theta_vec2.')';  % [Ny^2 x K]
fr = F2 * Imag(S.Theta2);
scale = S.rX ./ mean(fr);
F2 = F2 * scale;

% ---- Tile into F exactly like the original triple-loop (but safe & general)
% Original layout:
%   - Rows: for each i in 1:Ny, there are 2*Ny rows (first Ny for F1, next Ny for F2)
%   - Cols: first K columns for F1, next K columns for F2
F = zeros(2*Ny*Ny, 2*K);
for i = 1:Ny
    for j = 1:Ny
        % row indices in the TILED matrix
        baseRow = (i-1)*(2*Ny) + j;         % position within the 2*Ny block
        % column blocks
        F(baseRow,             1:K)       = F1((j-1)*Ny + i, :);
        F(baseRow + Ny,        K+(1:K))   = F2((j-1)*Ny + i, :);
    end
end

end  % main function

% -------------------------------------------------------------------------
function theta_map = resolveOrientationMap(inputArg, varName)
% Accepts a numeric array or a .mat file path (char/string)
if isnumeric(inputArg)
    theta_map = inputArg;
    return;
end
if ~(ischar(inputArg) || isstring(inputArg))
    error('Invalid orientation map input. Provide a numeric matrix or a .mat file path.');
end
filePath = char(inputArg);
if ~isfile(filePath)
    % allow loading from MATLAB path without extension
    [~,~,ext] = fileparts(filePath);
    if isempty(ext)
        if isfile([filePath '.mat'])
            filePath = [filePath '.mat'];
        end
    end
end
if ~isfile(filePath)
    error('Could not find file: %s', filePath);
end

S = load(filePath, varName);
if ~isfield(S, varName)
    error('Variable "%s" not found in %s. Use Name-Value pair ''VarName?'' to specify.', varName, filePath);
end
theta_map = S.(varName);
end
