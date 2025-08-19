function [Wrr,Wrf,outdegreeEE,outdegreeIE,outdegreeEI,outdegreeII,outdegreeEX,outdegreeIX,seqindX,seqindE,seqindI,presynIDE,presynIDI,presynIDX] = ...
    generate_weights_broad_indegree(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,P_ts,TS_th,b1,varargin)
% GENERATE_WEIGHTS_BROAD_INDEGREE  Connectivity with broad (Gamma-like) indegree distribution.
%
% Syntax
%   [Wrr,Wrf,outdegreeEE,outdegreeIE,outdegreeEI,outdegreeII,outdegreeEX,outdegreeIX,seqindX,seqindE,seqindI,presynIDE,presynIDI,presynIDX] = ...
%       generate_weights_broad_indegree(Ne,Ni,Nx,sigmaRX,sigmaRR,Prr,Prx,P_ts,TS_th,b1, Name,Value, ...)
%
% Description
%   Version with per-postsynaptic indegree drawn from a Gamma-like distribution.
%
% Name-Value (optional)
%   'Seed','ThetaV1','ThetaE','ThetaV1Var','ThetaEVar' : see other weight functions.
%
% Outputs
%   Wrr, Wrf, degree summaries, sequence indices, and presynaptic index lists per cell.
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
validateattributes(b1,      {'double','single'},{'scalar','>',0}, mfilename, 'b1');

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

presynIDE = cell(1,Ne+Ni);
presynIDI = cell(1,Ne+Ni);
presynIDX = cell(1,Ne+Ni);

[I0,I1] = resolveTuning(NV.ThetaV1,NV.ThetaV1Var,Nx, NV.ThetaE,NV.ThetaEVar,Ne);

% ----- Wrr: E->E with Gamma-like indegree per E target
Wee2 = zeros(Ne*(Kee_in+Kee_ts_in),2);
alpha = (Kee_in+Kee_ts_in)/b1;
Kee_list     = gen_gamma_integer_array(alpha, b1, Ne, size(Wee2,1)); % total E->E per E
Kee_list_ts  = floor(Kee_list * P_ts(2));
Kee_list_ind = Kee_list - Kee_list_ts;

counter = 0;
for j=1:Ne
    x_pre = ceil(j/Ney); y_pre = mod(j-1,Ney)+1;
    x_post = CircRandN(x_pre, betaee, 1, Nex, Kee_list_ind(j));
    y_post = CircRandN(y_pre, betaee, 1, Ney, Kee_list_ind(j));
    nBlock = Kee_list_ind(j) + Kee_list_ts(j);
    Wee2(counter+(1:nBlock),1) = j;
    Wee2(counter+(1:Kee_list_ind(j)),2) = sort((x_post-1)*Ney + y_post);
    TS = cos((I1(j)-I1)*2*pi);
    idx = find(TS>TS_th); if isempty(idx), idx = (1:Ne)'; end
    Wee2(counter+Kee_list_ind(j)+(1:Kee_list_ts(j)),2) = randsample(idx, Kee_list_ts(j), true);
    presynIDE{j} = Wee2(counter+(1:nBlock),2);
    counter = counter + nBlock;
end

% ----- Wrr: I->E
Wei2 = zeros(Ne*Kei_in,2);
alpha = Kei_in/b1;  % nominal mean alpha*b1
Kei_list = gen_gamma_integer_array(alpha, b1, Ne, size(Wei2,1));
counter=0;
for j=1:Ne
    x_pre=ceil(j/Ney)*Niy/Ney; y_pre=(mod(j-1,Ney)+1)*Niy/Ney;
    x_post=CircRandN(x_pre,betaie,1,Nix,Kei_list(j));
    y_post=CircRandN(y_pre,betaie,1,Niy,Kei_list(j));
    Wei2(counter+(1:Kei_list(j)),1)=j;
    Wei2(counter+(1:Kei_list(j)),2)=sort((x_post-1)*Niy+y_post+Ne);
    presynIDI{j}=Wei2(counter+(1:Kei_list(j)),2);
    counter=counter+Kei_list(j);
end

% ----- Wrr: E->I
Wie2 = zeros(Ni*Kie_in,2);
alpha = Kie_in/b1;
Kie_list = gen_gamma_integer_array(alpha, b1, Ni, size(Wie2,1));
counter=0;
for j=1:Ni
    x_pre=ceil(j/Niy)*Ney/Niy; y_pre=(mod(j-1,Niy)+1)*Ney/Niy;
    x_post=CircRandN(x_pre,betaei,1,Nex,Kie_list(j));
    y_post=CircRandN(y_pre,betaei,1,Ney,Kie_list(j));
    Wie2(counter+(1:Kie_list(j)),1)=Ne+j;
    Wie2(counter+(1:Kie_list(j)),2)=sort((x_post-1)*Ney+y_post);
    presynIDE{j+Ne}=Wie2(counter+(1:Kie_list(j)),2);
    counter=counter+Kie_list(j);
end

% ----- Wrr: I->I
Wii2 = zeros(Ni*Kii_in,2);
alpha = Kii_in/b1;
Kii_list = gen_gamma_integer_array(alpha, b1, Ni, size(Wii2,1));
counter=0;
for j=1:Ni
    x_pre=ceil(j/Niy); y_pre=(mod(j-1,Niy)+1);
    x_post=CircRandN(x_pre,betaii,1,Nix,Kii_list(j));
    y_post=CircRandN(y_pre,betaii,1,Niy,Kii_list(j));
    Wii2(counter+(1:Kii_list(j)),1)=Ne+j;
    Wii2(counter+(1:Kii_list(j)),2)=sort((x_post-1)*Niy+y_post+Ne);
    presynIDI{j+Ne}=Wii2(counter+(1:Kii_list(j)),2);
    counter=counter+Kii_list(j);
end

% Collect Wrr by postsynaptic target
jj=0; outdegreeEE=zeros(Ne,1); outdegreeIE=zeros(Ne,1);
for j=1:Ne
    outdegreeEE(j)=outdegreeEE(j)+sum(Wee2(:,2)==j);
    tempseq=Wee2(Wee2(:,2)==j,1);
    Wrr(jj+(1:numel(tempseq)))=tempseq; jj=jj+numel(tempseq);
    outdegreeIE(j)=outdegreeIE(j)+sum(Wie2(:,2)==j);
    tempseq=Wie2(Wie2(:,2)==j,1);
    Wrr(jj+(1:numel(tempseq)))=tempseq; jj=jj+numel(tempseq);
end
outdegreeE=outdegreeEE+outdegreeIE; seqindE=cumsum(outdegreeE); seqindE=[0;seqindE];

outdegreeEI=zeros(Ni,1); outdegreeII=zeros(Ni,1);
for j=1:Ni
    outdegreeEI(j)=outdegreeEI(j)+sum(Wei2(:,2)==j+Ne);
    tempseq=Wei2(Wei2(:,2)==j+Ne,1);
    Wrr(jj+(1:numel(tempseq)))=tempseq; jj=jj+numel(tempseq);
    outdegreeII(j)=outdegreeII(j)+sum(Wii2(:,2)==j+Ne);
    tempseq=Wii2(Wii2(:,2)==j+Ne,1);
    Wrr(jj+(1:numel(tempseq)))=tempseq; jj=jj+numel(tempseq);
end
outdegreeI=outdegreeEI+outdegreeII; seqindI=cumsum(outdegreeI); seqindI=[0;seqindI];

% ----- Wrf
Wex2 = zeros(Ne*(Kex_in+Kex_ts_in),2);
alpha = (Kex_in+Kex_ts_in)/b1;
Kex_list = gen_gamma_integer_array(alpha, b1, Ne, size(Wex2,1));
Kex_list_ts = floor(Kex_list*P_ts(1));
Kex_list_ind = Kex_list - Kex_list_ts;

counter=0;
for j=1:Ne
    x_pre=ceil(j/Ney)*Nxy/Ney; y_pre=(mod(j-1,Ney)+1)*Nxy/Ney;
    x_post=CircRandN(x_pre,betaex,1,Nxx,Kex_list_ind(j));
    y_post=CircRandN(y_pre,betaex,1,Nxy,Kex_list_ind(j));
    nBlock = Kex_list_ind(j)+Kex_list_ts(j);
    Wex2(counter+(1:nBlock),1)=j;
    Wex2(counter+(1:Kex_list_ind(j)),2)=sort((x_post-1)*Nxy+y_post);
    TS=cos((I1(j)-I0)*2*pi);
    idx=find(TS>TS_th); if isempty(idx), idx=(1:Nx)'; end
    Wex2(counter+Kex_list_ind(j)+(1:Kex_list_ts(j)),2)=randsample(idx,Kex_list_ts(j),true);
    presynIDX{j}=Wex2(counter+(1:nBlock),2);
    counter=counter+nBlock;
end

Wix2 = zeros(Ni*Kix_in,2);
alpha = Kix_in/b1;
Kix_list = gen_gamma_integer_array(alpha, b1, Ni, size(Wix2,1));
counter=0;
for j=1:Ni
    x_pre=ceil(j/Niy)*Nxy/Niy; y_pre=(mod(j-1,Niy)+1)*Nxy/Niy;
    x_post=CircRandN(x_pre,betaix,1,Nxx,Kix_list(j));
    y_post=CircRandN(y_pre,betaix,1,Nxy,Kix_list(j));
    Wix2(counter+(1:Kix_list(j)),1)=j+Ne;
    Wix2(counter+(1:Kix_list(j)),2)=sort((x_post-1)*Nxy+y_post);
    presynIDX{j+Ne}=Wix2(counter+(1:Kix_list(j)),2);
    counter=counter+Kix_list(j);
end

jj=0; outdegreeEX=zeros(Nx,1); outdegreeIX=zeros(Nx,1);
for j=1:Nx
    outdegreeEX(j)=outdegreeEX(j)+sum(Wex2(:,2)==j);
    tempseq=Wex2(Wex2(:,2)==j,1);
    Wrf(jj+(1:numel(tempseq)))=tempseq; jj=jj+numel(tempseq);
    outdegreeIX(j)=outdegreeIX(j)+sum(Wix2(:,2)==j);
    tempseq=Wix2(Wix2(:,2)==j,1);
    Wrf(jj+(1:numel(tempseq)))=tempseq; jj=jj+numel(tempseq);
end
outdegreeX=outdegreeEX+outdegreeIX; seqindX=cumsum(outdegreeX); seqindX=[0;seqindX];

end % function

% --------------------------- Local helpers --------------------------------
function validateCounts(Ne,Ni,Nx); for v={'Ne','Ni','Nx'}; name=v{1}; val=eval(name);
if ~(isnumeric(val)&&isscalar(val)&&val>0&&mod(val,1)==0), error('Argument %s must be a positive integer.',name); end
end; end

function [Nex,Ney,Nix,Niy,Nxx,Nxy] = inferGridDims(Ne,Ni,Nx)
Nex=sqrt(Ne/2); Ney=Nex*2; Nix=sqrt(Ni/2); Niy=Nix*2; Nxx=sqrt(Nx/2); Nxy=Nxx*2;
if any(abs([Nex,Nix,Nxx]-round([Nex,Nix,Nxx]))>1e-9)
    error('Ne, Ni, Nx must satisfy Ne=2*n^2 etc. Got Ne=%d Ni=%d Nx=%d',Ne,Ni,Nx);
end
Nex=round(Nex); Nix=round(Nix); Nxx=round(Nxx);
end

