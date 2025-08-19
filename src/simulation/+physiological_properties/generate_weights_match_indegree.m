function [Wrr,Wrf,outdegreeEE,outdegreeIE,outdegreeEI,outdegreeII,outdegreeEX,outdegreeIX,seqindX,seqindE,seqindI] = ...
    generate_weights_match_indegree(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,P_ts,TS_th,varargin)
% GENERATE_WEIGHTS_MATCH_INDEGREE  Connectivity with matched indegree per postsynaptic unit.
%
% Syntax
%   [Wrr,Wrf,outdegreeEE,outdegreeIE,outdegreeEI,outdegreeII,outdegreeEX,outdegreeIX,seqindX,seqindE,seqindI] = ...
%       generate_weights_match_indegree(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,P_ts,TS_th, Name,Value, ...)
%
% Description
%   Similar to GENERATE_WEIGHTS_DEFAULT, but constructs auxiliary tables (outdegree
%   and sequence indices) to match indegree constraints while preserving spatial and
%   tuning structure. Refactored to remove hard-coded paths and improve robustness.
%
% Name-Value (optional)
%   'Seed'       : rng seed.
%   'ThetaV1'    : V1 orientation map (numeric or file path). See default version.
%   'ThetaE'     : E orientation vector/map (numeric or file path). See default.
%   'ThetaV1Var' : variable name for ThetaV1 file (default 'theta_map').
%   'ThetaEVar'  : variable name for ThetaE file  (default 'theta_mapE').
%
% Outputs
%   Wrr, Wrf : connectivity lists (int32).
%   outdegree??, seqind?? : degree summaries and sequence indices.
%
% Example
%   [Wrr,Wrf,~,~,~,~,~,~,seqX,seqE,seqI] = generate_weights_match_indegree(...,'Seed',1);
%
% Author: Refactored for clarity & robustness.
% -------------------------------------------------------------------------

p = inputParser; p.FunctionName = mfilename;
addParameter(p,'Seed',[],@(x)isnumeric(x)&&isscalar(x) || isempty(x));
addParameter(p,'ThetaV1',[],@(x)isnumeric(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p,'ThetaE',[],@(x)isnumeric(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p,'ThetaV1Var','theta_map',@(s)ischar(s)||isstring(s));
addParameter(p,'ThetaEVar','theta_mapE',@(s)ischar(s)||isstring(s));
parse(p,varargin{:}); NV = p.Results;
if ~isempty(NV.Seed), rng(NV.Seed,'twister'); end

validateCounts(Ne,Ni,Nx);
validateattributes(sigmaRX, {'double','single'},{'size',[2,1]}, mfilename, 'sigmaRX');
validateattributes(sigmaRR, {'double','single'},{'size',[2,2]}, mfilename, 'sigmaRR');
validateattributes(Prr,     {'double','single'},{'size',[2,2]}, mfilename, 'Prr');
validateattributes(Prx,     {'double','single'},{'size',[2,1]}, mfilename, 'Prx');
validateattributes(P_ts,    {'double','single'},{'size',[2,1],'>=',0,'<=',1}, mfilename, 'P_ts');
validateattributes(TS_th,   {'double','single'},{'scalar'}, mfilename, 'TS_th');

sigmaeX=sigmaRX(1); sigmaiX=sigmaRX(2);
sigmaee=sigmaRR(1,1); sigmaei=sigmaRR(1,2);
sigmaie=sigmaRR(2,1); sigmaii=sigmaRR(2,2);

pee0=Prr(1,1); pei0=Prr(1,2); pie0=Prr(2,1); pii0=Prr(2,2);
pex0=Prx(1);   pix0=Prx(2);

[Nex,Ney,Nix,Niy,Nxx,Nxy] = inferGridDims(Ne,Ni,Nx);

betaee=sigmaee*Nex; betaei=sigmaei*Nex; betaie=sigmaie*Nix; betaii=sigmaii*Nix;
betaex=sigmaeX*Nex; betaix=sigmaiX*Nix;

Kee_in   = ceil(pee0*Ne*(1-P_ts(1)));
Kee_ts_in= floor(pee0*Ne*(P_ts(1)));
Kei_in   = ceil(pei0*Ni);
Kie_in   = ceil(pie0*Ne);
Kii_in   = ceil(pii0*Ni);
Kex_in   = ceil(pex0*Nx*(1-P_ts(2)));
Kex_ts_in= floor(pex0*Nx*(P_ts(2)));
Kix_in   = ceil(pix0*Nx);

CircRandN=@(mu,sigma,minv,maxv,n)(mod(round(sigma*randn(n,1)+mu)-minv, maxv-minv+1)+minv);

Ke_in = Kee_in + Kee_ts_in + Kei_in;
Ki_in = Kie_in + Kii_in;
Wrr = zeros(Ke_in*Ne + Ki_in*Ni, 1, 'int32');
Wrf = zeros((Kex_in+Kex_ts_in)*Ne + Kix_in*Ni, 1, 'int32');

[I0,I1] = resolveTuning(NV.ThetaV1,NV.ThetaV1Var,Nx, NV.ThetaE,NV.ThetaEVar,Ne);

% ----- Wrr
Wee2 = zeros(Ne*(Kee_in+Kee_ts_in),2);
for j=1:Ne
    x_pre=ceil(j/Ney); y_pre=mod(j-1,Ney)+1;
    x_post=CircRandN(x_pre,betaee,1,Nex,Kee_in);
    y_post=CircRandN(y_pre,betaee,1,Ney,Kee_in);
    idx0=(j-1)*(Kee_in+Kee_ts_in)+1;
    Wee2(idx0:idx0+Kee_in-1,1)=j;
    Wee2(idx0:idx0+Kee_in-1,2)=sort((x_post-1)*Ney+y_post);
    TS=cos((I1(j)-I1)*2*pi);
    idx=find(TS>TS_th); if isempty(idx), idx=1:Ne; end
    Wee2(idx0+Kee_in:idx0+Kee_in+Kee_ts_in-1,1)=j;
    Wee2(idx0+Kee_in:idx0+Kee_in+Kee_ts_in-1,2)=randsample(idx,Kee_ts_in,true);
end

Wei2 = zeros(Ne*Kei_in,2);
for j=1:Ne
    x_pre=ceil(j/Ney)*Niy/Ney; y_pre=(mod(j-1,Ney)+1)*Niy/Ney;
    x_post=CircRandN(x_pre,betaie,1,Nix,Kei_in);
    y_post=CircRandN(y_pre,betaie,1,Niy,Kei_in);
    idx0=(j-1)*Kei_in+1;
    Wei2(idx0:idx0+Kei_in-1,1)=j;
    Wei2(idx0:idx0+Kei_in-1,2)=sort((x_post-1)*Niy+y_post+Ne);
end

Wie2 = zeros(Ni*Kie_in,2);
for j=1:Ni
    x_pre=ceil(j/Niy)*Ney/Niy; y_pre=(mod(j-1,Niy)+1)*Ney/Niy;
    x_post=CircRandN(x_pre,betaei,1,Nex,Kie_in);
    y_post=CircRandN(y_pre,betaei,1,Ney,Kie_in);
    idx0=(j-1)*Kie_in+1;
    Wie2(idx0:idx0+Kie_in-1,1)=Ne+j;
    Wie2(idx0:idx0+Kie_in-1,2)=sort((x_post-1)*Ney+y_post);
end

Wii2 = zeros(Ni*Kii_in,2);
for j=1:Ni
    x_pre=ceil(j/Niy); y_pre=(mod(j-1,Niy)+1);
    x_post=CircRandN(x_pre,betaii,1,Nix,Kii_in);
    y_post=CircRandN(y_pre,betaii,1,Niy,Kii_in);
    idx0=(j-1)*Kii_in+1;
    Wii2(idx0:idx0+Kii_in-1,1)=Ne+j;
    Wii2(idx0:idx0+Kii_in-1,2)=sort((x_post-1)*Niy+y_post+Ne);
end

jj=0;
outdegreeEE=zeros(Ne,1); outdegreeIE=zeros(Ne,1);
for j=1:Ne
    outdegreeEE(j)=outdegreeEE(j)+sum(Wee2(:,2)==j);
    tempseq=Wee2(Wee2(:,2)==j,1);
    Wrr(jj+1:jj+numel(tempseq))=tempseq; jj=jj+numel(tempseq);
    outdegreeIE(j)=outdegreeIE(j)+sum(Wie2(:,2)==j);
    tempseq=Wie2(Wie2(:,2)==j,1);
    Wrr(jj+1:jj+numel(tempseq))=tempseq; jj=jj+numel(tempseq);
end
outdegreeE=outdegreeEE+outdegreeIE;
seqindE=cumsum(outdegreeE); seqindE=[0;seqindE];

outdegreeEI=zeros(Ni,1); outdegreeII=zeros(Ni,1);
for j=1:Ni
    outdegreeEI(j)=outdegreeEI(j)+sum(Wei2(:,2)==j+Ne);
    tempseq=Wei2(Wei2(:,2)==j+Ne,1);
    Wrr(jj+1:jj+numel(tempseq))=tempseq; jj=jj+numel(tempseq);
    outdegreeII(j)=outdegreeII(j)+sum(Wii2(:,2)==j+Ne);
    tempseq=Wii2(Wii2(:,2)==j+Ne,1);
    Wrr(jj+1:jj+numel(tempseq))=tempseq; jj=jj+numel(tempseq);
end
outdegreeI=outdegreeEI+outdegreeII;
seqindI=cumsum(outdegreeI); seqindI=[0;seqindI];

% ----- Wrf
Wex2 = zeros(Ne*(Kex_in+Kex_ts_in),2);
for j=1:Ne
    x_pre=ceil(j/Ney)*Nxy/Ney; y_pre=(mod(j-1,Ney)+1)*Nxy/Ney;
    x_post=CircRandN(x_pre,betaex,1,Nxx,Kex_in);
    y_post=CircRandN(y_pre,betaex,1,Nxy,Kex_in);
    idx0=(j-1)*(Kex_in+Kex_ts_in)+1;
    Wex2(idx0:idx0+Kex_in-1,1)=j;
    Wex2(idx0:idx0+Kex_in-1,2)=sort((x_post-1)*Nxy+y_post);
    TS=cos((I1(j)-I0)*2*pi);
    idx=find(TS>TS_th); if isempty(idx), idx=1:Nx; end
    Wex2(idx0+Kex_in:idx0+Kex_in+Kex_ts_in-1,1)=j;
    Wex2(idx0+Kex_in:idx0+Kex_in+Kex_ts_in-1,2)=randsample(idx,Kex_ts_in,true);
end

Wix2 = zeros(Ni*Kix_in,2);
for j=1:Ni
    x_pre=ceil(j/Niy)*Nxy/Niy; y_pre=(mod(j-1,Niy)+1)*Nxy/Niy;
    x_post=CircRandN(x_pre,betaix,1,Nxx,Kix_in);
    y_post=CircRandN(y_pre,betaix,1,Nxy,Kix_in);
    idx0=(j-1)*Kix_in+1;
    Wix2(idx0:idx0+Kix_in-1,1)=j+Ne;
    Wix2(idx0:idx0+Kix_in-1,2)=sort((x_post-1)*Nxy+y_post);
end

jj=0; outdegreeEX=zeros(Nx,1); outdegreeIX=zeros(Nx,1);
for j=1:Nx
    outdegreeEX(j)=outdegreeEX(j)+sum(Wex2(:,2)==j);
    tempseq=Wex2(Wex2(:,2)==j,1);
    Wrf(jj+1:jj+numel(tempseq))=tempseq; jj=jj+numel(tempseq);
    outdegreeIX(j)=outdegreeIX(j)+sum(Wix2(:,2)==j);
    tempseq=Wix2(Wix2(:,2)==j,1);
    Wrf(jj+1:jj+numel(tempseq))=tempseq; jj=jj+numel(tempseq);
end
outdegreeX=outdegreeEX+outdegreeIX;
seqindX=cumsum(outdegreeX); seqindX=[0;seqindX];

end  % function

% Local helpers reused from default version --------------------------------
function validateCounts(Ne,Ni,Nx)
for v = {'Ne','Ni','Nx'}
    name = v{1}; val = eval(name);
    if ~(isnumeric(val) && isscalar(val) && val>0 && mod(val,1)==0)
        error('Argument %s must be a positive integer.', name);
    end
end
end

function [Nex,Ney,Nix,Niy,Nxx,Nxy] = inferGridDims(Ne,Ni,Nx)
Nex = sqrt(Ne/2); Ney = Nex*2;
Nix = sqrt(Ni/2); Niy = Nix*2;
Nxx = sqrt(Nx/2); Nxy = Nxx*2;
if any(abs([Nex,Nix,Nxx]-round([Nex,Nix,Nxx])) > 1e-9)
    error(['Ne, Ni, Nx must satisfy Ne=2*n^2, Ni=2*n^2, Nx=2*n^2 for integer n. ' ...
           'Got Ne=%d, Ni=%d, Nx=%d.'], Ne, Ni, Nx);
end
Nex = round(Nex); Nix = round(Nix); Nxx = round(Nxx);
end

function [I0, I1] = resolveTuning(thetaV1, varV1, Nx, thetaE, varE, Ne)
if isempty(thetaV1) || isempty(thetaE)
    I0 = zeros(Nx,1); I1 = zeros(Ne,1); return;
end
if isnumeric(thetaV1), mapV1 = thetaV1; else
    S = load(char(thetaV1), varV1);
    if ~isfield(S,varV1), error('Variable "%s" missing in %s.', varV1, char(thetaV1)); end
    mapV1 = S.(varV1);
end
if numel(mapV1)~=Nx
    [r,c]=size(mapV1);
    if r*c~=Nx, error('ThetaV1 size mismatch (expected %d).', Nx); end
end
I0 = reshape(mapV1.',[],1);
if isnumeric(thetaE), mapE = thetaE; else
    S = load(char(thetaE), varE);
    if ~isfield(S,varE), error('Variable "%s" missing in %s.', varE, char(thetaE)); end
    mapE = S.(varE);
end
if numel(mapE)~=Ne
    [r,c]=size(mapE);
    if r*c~=Ne, error('ThetaE size mismatch (expected %d).', Ne); end
    mapE = reshape(mapE.',[],1);
else
    mapE = reshape(mapE,[],1);
end
I1 = mapE;
end
