function [Wrr, Wrf] = generate_weights_default(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,P_ts,TS_th,varargin)
% GENERATE_WEIGHTS_DEFAULT  Recurrent & feedforward connectivity (spatial + tuning).
%
% Syntax
%   [Wrr, Wrf] = generate_weights_default(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,P_ts,TS_th, Name,Value, ...)
%
% Description
%   Constructs sparse connectivity lists for a 2D E/I sheet receiving input X,
%   with connection probabilities shaped by spatial spread (sigmaRX, sigmaRR) and
%   optional tuning dependence. This refactor removes hard-coded paths and adds
%   robust input validation and reproducibility controls.
%
% Required Inputs
%   Ne, Ni, Nx : Positive integers. Counts of Excitatory, Inhibitory, and Input units.
%   sigmaRX    : [2x1] [sigma_eX; sigma_iX] spatial widths (in units of grid cells).
%   sigmaRR    : [2x2] [sigma_ee sigma_ei; sigma_ie sigma_ii] spatial widths.
%   Prr        : [2x2] connection probabilities within recurrent layer.
%   Prx        : [2x1] connection probabilities from input X: [p_eX; p_iX].
%   P_ts       : [2x1] fraction of connections that are tuning-dependent:
%                   P_ts(1): E->E tuning-dependent fraction
%                   P_ts(2): X->E tuning-dependent fraction
%   TS_th      : Scalar threshold for tuning similarity (cosine on angles*2*pi).
%
% Name-Value Pairs (optional)
%   'Seed'         : Numeric scalar for rng(Seed,'twister'). Default: [] (no change).
%   'ThetaV1'      : Either a numeric [Nx x Nx] orientation map or a file path
%                    (string/char) to a .mat containing variable 'theta_map' (or
%                    specify name via 'ThetaV1Var'). Used for input X.
%   'ThetaE'       : Either a numeric [Ne x 1] orientation vector for E units,
%                    or a file path to a .mat with variable 'theta_mapE' (or
%                    specify via 'ThetaEVar'). If a 2D [Ny x Nx] map is provided,
%                    it will be reshaped column-wise to [Ne x 1] (with Ne = Ny*Nx).
%   'ThetaV1Var'   : Variable name in ThetaV1 file. Default 'theta_map'.
%   'ThetaEVar'    : Variable name in ThetaE file.  Default 'theta_mapE'.
%
% Outputs
%   Wrr : int32 vector of presynaptic indices for recurrent targets, concatenated
%         by postsynaptic cell (EE/EI then IE/II as in original design).
%   Wrf : int32 vector of presynaptic indices for feedforward connections (X->E/I).
%
% Example
%   [Wrr, Wrf] = generate_weights_default(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,[.2;.3],0.2, ...
%                  'Seed',1,'ThetaV1','theta_map_V1_two_rec.mat','ThetaE','theta_map_MTE_two_rec.mat');
%
% Author: Refactored for clarity & robustness.
% -------------------------------------------------------------------------

% --- Parse Name-Value pairs
p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'Seed',[],@(x)isnumeric(x)&&isscalar(x) || isempty(x));
addParameter(p,'ThetaV1',[],@(x)isnumeric(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p,'ThetaE',[],@(x)isnumeric(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p,'ThetaV1Var','theta_map',@(s)ischar(s)||isstring(s));
addParameter(p,'ThetaEVar','theta_mapE',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
NV = p.Results;

if ~isempty(NV.Seed)
    rng(NV.Seed,'twister');
end

% --- Validate core numeric inputs
validateCounts(Ne,Ni,Nx);
validateattributes(sigmaRX, {'double','single'},{'size',[2,1]}, mfilename, 'sigmaRX');
validateattributes(sigmaRR, {'double','single'},{'size',[2,2]}, mfilename, 'sigmaRR');
validateattributes(Prr,     {'double','single'},{'size',[2,2]}, mfilename, 'Prr');
validateattributes(Prx,     {'double','single'},{'size',[2,1]}, mfilename, 'Prx');
validateattributes(P_ts,    {'double','single'},{'size',[2,1],'>=',0,'<=',1}, mfilename, 'P_ts');
validateattributes(TS_th,   {'double','single'},{'scalar'}, mfilename, 'TS_th');

% --- Unpack spreads
sigmaeX = sigmaRX(1); sigmaiX = sigmaRX(2);
sigmaee = sigmaRR(1,1); sigmaei = sigmaRR(1,2);
sigmaie = sigmaRR(2,1); sigmaii = sigmaRR(2,2);

pee0 = Prr(1,1); pei0 = Prr(1,2);
pie0 = Prr(2,1); pii0 = Prr(2,2);
pex0 = Prx(1);   pix0 = Prx(2);

% --- Grid sizes (require square tilings Ne=2*Nex^2 etc., as in original code)
[Nex,Ney,Nix,Niy,Nxx,Nxy] = inferGridDims(Ne,Ni,Nx);

betaee = sigmaee*Nex;
betaei = sigmaei*Nex;
betaie = sigmaie*Nix;
betaii = sigmaii*Nix;
betaex = sigmaeX*Nex;
betaix = sigmaiX*Nix;

% --- Connection counts
Kee    = ceil(pee0*Ne*(1-P_ts(1)));
Kee_ts = ceil(pee0*Ne*(P_ts(1)));
Kei    = ceil(pei0*Ne);
Kie    = ceil(pie0*Ni);
Kii    = ceil(pii0*Ni);
Kex    = ceil(pex0*Ne*(1-P_ts(2)));
Kex_ts = ceil(pex0*Ne*(P_ts(2)));
Kix    = ceil(pix0*Ni);

CircRandN = @(mu,sigma,minv,maxv,n)(mod(round(sigma*randn(n,1)+mu)-minv, maxv-minv+1)+minv);

Ke = Kee+Kee_ts+Kie;
Ki = Kei+Kii;
Kx = Kex+Kex_ts+Kix;

Wrr = zeros(Ke*Ne + Ki*Ni, 1, 'int32');
Wrf = zeros(Kx*Nx, 1, 'int32');

% --- Resolve tuning maps (optional; only needed for tuning-dependent parts)
[I0, I1] = resolveTuning(NV.ThetaV1,NV.ThetaV1Var,Nx, NV.ThetaE,NV.ThetaEVar,Ne);

% =================== Build Wrr ===================
for j = 1:Ne
    % E pre, E post (tuning-independent)
    x_pre = ceil(j/Ney);  y_pre = mod(j-1,Ney)+1;
    x_post = CircRandN(x_pre, betaee, 1, Nex, Kee);
    y_post = CircRandN(y_pre, betaee, 1, Ney, Kee);
    Wrr((1+(j-1)*Ke):(Kee+(j-1)*Ke)) = sort((x_post-1)*Ney + y_post);

    % E pre, E post (tuning-dependent)
    if Kee_ts>0
        TS = cos((I1(j) - I1)*2*pi);
        idx = find(TS > TS_th);
        if isempty(idx), idx = 1:Ne; end  % graceful fallback
        Wrr((Kee+1+(j-1)*Ke):(Kee+Kee_ts+(j-1)*Ke)) = randsample(idx, Kee_ts, true);
    end

    % E pre, I post
    x_pre = ceil(j/Ney)*Niy/Ney;  y_pre = (mod(j-1,Ney)+1)*Niy/Ney;
    x_post = CircRandN(x_pre, betaie, 1, Nix, Kie);
    y_post = CircRandN(y_pre, betaie, 1, Niy, Kie);
    Wrr((Kee+Kee_ts+1+(j-1)*Ke):(j*Ke)) = sort((x_post-1)*Niy + y_post + Ne);
end

offset = Ne*Ke;
for j = 1:Ni
    % I pre, E post
    x_pre = ceil(j/Niy)*Ney/Niy;  y_pre = (mod(j-1,Niy)+1)*Ney/Niy;
    x_post = CircRandN(x_pre, betaei, 1, Nex, Kei);
    y_post = CircRandN(y_pre, betaei, 1, Ney, Kei);
    Wrr((offset+1+(j-1)*Ki):(offset+Kei+(j-1)*Ki)) = sort((x_post-1)*Ney + y_post);

    % I pre, I post
    x_pre = ceil(j/Niy);  y_pre = (mod(j-1,Niy)+1);
    x_post = CircRandN(x_pre, betaii, 1, Nix, Kii);
    y_post = CircRandN(y_pre, betaii, 1, Niy, Kii);
    Wrr((offset+Kei+1+(j-1)*Ki):(offset+Kei+Kii+(j-1)*Ki)) = sort((x_post-1)*Niy + y_post + Ne);
end

% =================== Build Wrf ===================
for j = 1:Nx
    % X pre, E post (tuning-independent)
    x_pre = ceil(j/Nxy)*Ney/Nxy;  y_pre = (mod(j-1,Nxy)+1)*Ney/Nxy;
    x_post = CircRandN(x_pre, betaex, 1, Nex, Kex);
    y_post = CircRandN(y_pre, betaex, 1, Ney, Kex);
    Wrf((Kx*(j-1)+1):(Kx*(j-1)+Kex)) = sort((x_post-1)*Ney + y_post);

    % X pre, E post (tuning-dependent)
    if Kex_ts>0
        TS = cos((I0(j) - I1)*2*pi);
        idx = find(TS > TS_th);
        if isempty(idx), idx = 1:Ne; end
        Wrf((Kex+1+(j-1)*Kx):(Kex+Kex_ts+(j-1)*Kx)) = randsample(idx, Kex_ts, true);
    end

    % X pre, I post
    x_pre = ceil(j/Nxy)*Niy/Nxy;  y_pre = (mod(j-1,Nxy)+1)*Niy/Nxy;
    x_post = CircRandN(x_pre, betaix, 1, Nix, Kix);
    y_post = CircRandN(y_pre, betaix, 1, Niy, Kix);
    Wrf((Kex+Kex_ts+1+(j-1)*Kx):(j*Kx)) = sort((x_post-1)*Niy + y_post + Ne);
end

end % function

% -------------------------------------------------------------------------
function validateCounts(Ne,Ni,Nx)
% basic sanity checks
for v = {'Ne','Ni','Nx'}
    name = v{1};
    val = eval(name);
    if ~(isnumeric(val) && isscalar(val) && val>0 && mod(val,1)==0)
        error('Argument %s must be a positive integer.', name);
    end
end
end

% -------------------------------------------------------------------------
function [Nex,Ney,Nix,Niy,Nxx,Nxy] = inferGridDims(Ne,Ni,Nx)
% Recover sheet tilings used by the original code: Ne = 2*(Nex)^2, etc.
Nex = sqrt(Ne/2); Ney = Nex*2;
Nix = sqrt(Ni/2); Niy = Nix*2;
Nxx = sqrt(Nx/2); Nxy = Nxx*2;
if any(abs([Nex,Nix,Nxx]-round([Nex,Nix,Nxx])) > 1e-9)
    error(['Ne, Ni, Nx must satisfy Ne=2*n^2, Ni=2*n^2, Nx=2*n^2 for integer n. ' ...
           'Got Ne=%d, Ni=%d, Nx=%d.'], Ne, Ni, Nx);
end
Nex = round(Nex); Nix = round(Nix); Nxx = round(Nxx);
end

% -------------------------------------------------------------------------
function [I0, I1] = resolveTuning(thetaV1, varV1, Nx, thetaE, varE, Ne)
% Resolve/validate tuning vectors/maps used in tuning-dependent sampling
% I0(j): input orientation for X unit j (length Nx)
% I1(j): orientation for E unit j       (length Ne)
if isempty(thetaV1) || isempty(thetaE)
    % Allow operation without tuning (then tuning-dependent parts fall back gracefully)
    I0 = zeros(Nx,1);
    I1 = zeros(Ne,1);
    return;
end

% V1 (input) map
if isnumeric(thetaV1)
    mapV1 = thetaV1;
else
    S = load(char(thetaV1), varV1);
    if ~isfield(S, varV1), error('Variable "%s" not found in file %s.', varV1, char(thetaV1)); end
    mapV1 = S.(varV1);
end
if ~isnumeric(mapV1) || numel(mapV1)~=Nx
    % allow 2D map reshape
    if numel(mapV1)~=Nx
        [r,c] = size(mapV1);
        if r*c ~= Nx
            error('ThetaV1 size mismatch: expected Nx=%d elements, got %d.', Nx, numel(mapV1));
        end
    end
end
I0 = reshape(mapV1.', [], 1); % column vector of length Nx

% E map
if isnumeric(thetaE)
    mapE = thetaE;
else
    S = load(char(thetaE), varE);
    if ~isfield(S, varE), error('Variable "%s" not found in file %s.', varE, char(thetaE)); end
    mapE = S.(varE);
end
if numel(mapE) ~= Ne
    [r,c] = size(mapE);
    if r*c ~= Ne
        error('ThetaE size mismatch: expected Ne=%d elements, got %d.', Ne, numel(mapE));
    end
    mapE = reshape(mapE.', [], 1);
else
    mapE = reshape(mapE, [], 1);
end

I1 = mapE;
end
