function out = gen_V1_FiringRates(theta1, theta2, c1, c2, varargin)
%GENERATEV1FIRINGRATES  Build V1 firing-rate fields from filters + (Gabor or spontaneous) inputs.
%
% Syntax
%   out = gen_V1_FiringRates(theta1, theta2, c1, c2)
%   out = gen_V1_FiringRates(..., 'Image1','gabor'|'spont', 'Image2','gabor'|'spont', ...
%                               'Sigma',0.2, 'Lambda',0.6, 'FilterFile','V1filterRecSig0d2Lam0d6.mat', ...
%                               'DataDir',pwd, 'OutDir',pwd, 'Save',true, 'BaseName','V1fr')
%
% Inputs
%   theta1, theta2 : orientation (radians) for left/right images (used if 'gabor')
%   c1, c2         : contrast scalar multipliers (used if 'gabor')
%
% Name-Value Pairs
%   'Image1','Image2' : 'gabor' or 'spont' (default 'gabor' for both)
%   'Sigma'           : Gabor Gaussian sigma (default 0.2)
%   'Lambda'          : Gabor wavelength (default 0.6)
%   'FilterFile'      : .mat file containing variable F (default 'V1filterRecSig0d2Lam0d6.mat')
%   'DataDir'         : folder to load FilterFile from (default: if figureConstants exists, C.dir.data; else pwd)
%   'OutDir'          : folder to save outputs (default: if figureConstants exists, C.dir.result; else pwd)
%   'Save'            : logical, save *_Rec_1/2/3 .mat files (default true)
%   'BaseName'        : filename prefix (default 'V1fr')
%
% Output (struct)
%   out.fr_both   : 10000x1 vector (left=Image1, right=Image2)
%   out.fr_left   : 10000x1 vector (left=Image1, right=5 Hz)
%   out.fr_right  : 10000x1 vector (left=5 Hz,    right=Image2)
%   out.meta      : parameters and paths used
%
% Notes
%   - Preserves original row/col packing via fr = reshape(FR', [], 1).
%   - Uses the same F → F1/F2 mapping as your nested loops, but vectorized.
%   - If 'spont', the corresponding 50x50 half is a homogeneous 5 Hz plane.
%
% Example
%   out = generateV1FiringRates(0.1*pi, 0.6*pi, 1.0, 0.8, 'Image1','gabor','Image2','spont');

  % -------------------- Parse inputs --------------------
  p = inputParser;
  addRequired(p, 'theta1', @(x) isscalar(x) && isnumeric(x));
  addRequired(p, 'theta2', @(x) isscalar(x) && isnumeric(x));
  addRequired(p, 'c1',     @(x) isscalar(x) && isnumeric(x));
  addRequired(p, 'c2',     @(x) isscalar(x) && isnumeric(x));
  addParameter(p, 'Image1', 'gabor', @(s) any(validatestring(s,{'gabor','spont'})));
  addParameter(p, 'Image2', 'gabor', @(s) any(validatestring(s,{'gabor','spont'})));
  addParameter(p, 'Sigma',  0.2,     @(x) isscalar(x) && x>0);
  addParameter(p, 'Lambda', 0.6,     @(x) isscalar(x) && x>0);
  addParameter(p, 'FilterFile', 'V1filterRecSig0d2Lam0d6.mat', @(s) ischar(s) || isstring(s));

  % Defaults for DataDir/OutDir from figureConstants if available
  thisFile = mfilename('fullpath');
  repoRoot = fileparts(fileparts(fileparts(fileparts(thisFile))));
  defaultOutDir = fullfile(repoRoot,'data');
  defaultDataDir = fullfile(repoRoot,'data');
  addParameter(p, 'DataDir', defaultDataDir, @(s) ischar(s) || isstring(s));
  addParameter(p, 'OutDir',  defaultOutDir,  @(s) ischar(s) || isstring(s));
  addParameter(p, 'Save',    true,           @(b) isscalar(b) && islogical(b));
  addParameter(p, 'BaseName','V1fr',         @(s) ischar(s) || isstring(s));
  parse(p, theta1, theta2, c1, c2, varargin{:});
  P = p.Results;

  % -------------------- Constants --------------------
  SPONT_RATE_HZ = 0.005;  % homogeneous spontaneous rate
  NX_HALF       = 50;     % 50x50 for each half (left/right)
  NPIX_HALF     = NX_HALF * NX_HALF;   % 2500
  NPIX_GABOR    = 225;    % 15x15 gabor image
  COLS_PER_SIDE = 225;    % F columns per side (matches NPIX_GABOR)
  ROWS_PER_I    = 100;    % F has 100 rows per i-block (top/bottom halves)
  c_base = 0.5;

  % -------------------- Load filter bank --------------------
  assert(isfolder(P.DataDir), 'DataDir not found: %s', P.DataDir);
  filtPath = fullfile(P.DataDir, P.FilterFile);
  assert(isfile(filtPath), 'Filter file not found: %s', filtPath);
  S = load(filtPath);
  assert(isfield(S,'F'), 'Variable "F" not found in %s', P.FilterFile);
  F = S.F;

  % Basic size sanity checks to match original layout expectations
  [nRows, nCols] = size(F);
  assert(mod(nRows, ROWS_PER_I)==0, 'F row count (%d) not divisible by 100.', nRows);
  nI = nRows / ROWS_PER_I;                 % expected 50
  assert(nI==NX_HALF, 'Expected 50 i-blocks; got %d (rows=%d).', nI, nRows);
  assert(nCols >= 2*COLS_PER_SIDE, 'F needs at least %d columns; got %d.', 2*COLS_PER_SIDE, nCols);

  % -------------------- Vectorized F → F1/F2 mapping --------------------
  % Left half uses rows j=1:50 across each 100-row block and cols 1:225
  Fblock1 = reshape(F(:, 1:COLS_PER_SIDE), [ROWS_PER_I, nI, COLS_PER_SIDE]);   % [100,50,225]
  F1 = reshape(permute(Fblock1(1:50,:,:), [2 1 3]), [NX_HALF*NX_HALF, COLS_PER_SIDE]);  % [2500,225]

  % Right half uses rows j=51:100 and cols 226:450
  Fblock2 = reshape(F(:, COLS_PER_SIDE+1:2*COLS_PER_SIDE), [ROWS_PER_I, nI, COLS_PER_SIDE]); % [100,50,225]
  F2 = reshape(permute(Fblock2(51:100,:,:), [2 1 3]), [NX_HALF*NX_HALF, COLS_PER_SIDE]);     % [2500,225]

  % -------------------- Build inputs: Gabor or spontaneous --------------------
  % 15x15 grid (matches NPIX_GABOR) centered at 0, spacing ≈ 0.07 (as in your script)
  side = sqrt(NPIX_GABOR);
  assert(abs(side - round(side)) < 1e-12, 'Gabor pixel count must be a perfect square.');
  side = round(side);
  x = linspace(-0.49, 0.49, side);
  [X, Y] = meshgrid(x, x);
  X = X(:); Y = Y(:);

  % Gabor image generator (column vector length 225)
  function g = gaborImage(theta, sigma, lambda)
    g = exp(-(X.^2 + Y.^2) / (2*sigma^2)) .* cos((2*pi/lambda) * (X*cos(theta) + Y*sin(theta)));
    g = g(:);  % ensure column
  end

  % Left image vector (225x1) or spontaneous
  if strcmpi(P.Image1, 'gabor')
    Imag1 = P.c1 * gaborImage(P.theta1, P.Sigma, P.Lambda) / c_base;
    assert(numel(Imag1) == COLS_PER_SIDE);
    frLeft  = reshape(F1 * Imag1, NX_HALF, NX_HALF);    % 50x50
  else
    frLeft  = SPONT_RATE_HZ * ones(NX_HALF, NX_HALF);
  end

  % Right image vector (225x1) or spontaneous
  if strcmpi(P.Image2, 'gabor')
    Imag2 = P.c2 * gaborImage(P.theta2, P.Sigma, P.Lambda) / c_base;
    assert(numel(Imag2) == COLS_PER_SIDE);
    frRight = reshape(F2 * Imag2, NX_HALF, NX_HALF);    % 50x50
  else
    frRight = SPONT_RATE_HZ * ones(NX_HALF, NX_HALF);
  end

  % -------------------- Assemble 100x100 fields and vectorize like original --------------------
  % Both halves present
  FR_both          = zeros(50,100);
  FR_both(1:50,1:50)   = frLeft;
  FR_both(1:50,51:100) = frRight;
  fr_both_vec      = reshape(FR_both', [], 1);

  % Left only (right spontaneous)
  FR_left          = zeros(50,100);
  FR_left(1:50,1:50)   = frLeft;
  FR_left(1:50,51:100) = SPONT_RATE_HZ * ones(NX_HALF, NX_HALF);
  fr_left_vec      = reshape(FR_left', [], 1);

  % Right only (left spontaneous)
  FR_right         = zeros(50,100);
  FR_right(1:50,1:50)   = SPONT_RATE_HZ * ones(NX_HALF, NX_HALF);
  FR_right(1:50,51:100) = frRight;
  fr_right_vec     = reshape(FR_right', [], 1);

  % -------------------- Save (match original naming) --------------------
  if P.Save
    if ~isfolder(P.OutDir), mkdir(P.OutDir); end
    fr = fr_both_vec;  save(fullfile(P.OutDir, sprintf('%s_Rec_3.mat', P.BaseName)), 'fr', '-v7'); 
    fr = fr_left_vec;  save(fullfile(P.OutDir, sprintf('%s_Rec_1.mat', P.BaseName)), 'fr', '-v7');  
    fr = fr_right_vec; save(fullfile(P.OutDir, sprintf('%s_Rec_2.mat', P.BaseName)), 'fr', '-v7'); 
  end

  % -------------------- Return struct --------------------
  out = struct();
  out.fr_both  = fr_both_vec;
  out.fr_left  = fr_left_vec;
  out.fr_right = fr_right_vec;
  out.meta = struct( ...
      'theta1',P.theta1,'theta2',P.theta2,'c1',P.c1,'c2',P.c2, ...
      'Image1',P.Image1,'Image2',P.Image2,'Sigma',P.Sigma,'Lambda',P.Lambda, ...
      'FilterFile',P.FilterFile,'DataDir',P.DataDir,'OutDir',P.OutDir,'Saved',P.Save);

end
