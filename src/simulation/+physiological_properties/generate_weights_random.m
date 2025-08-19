function [Wrr,Wrf] = generate_weights_random(Ne,Ni,Nx,Prr,Prx,varargin)
% GENERATE_WEIGHTS_RANDOM  Random (Erdos–Renyi) connectivity without spatial/tuning structure.
%
% Syntax
%   [Wrr, Wrf] = generate_weights_random(Ne,Ni,Nx,Prr,Prx, Name,Value, ...)
%
% Description
%   Constructs random connectivity lists according to global probabilities in Prr and Prx.
%   No tuning or spatial structure is used and no external files are loaded.
%
% Inputs
%   Ne, Ni, Nx : positive integers.
%   Prr : [2x2] probabilities within recurrent layer.
%   Prx : [2x1] probabilities from input X.
%
% Name-Value (optional)
%   'Seed' : rng seed for reproducibility.
%
% Outputs
%   Wrr, Wrf : int32 connectivity lists.
%
% Example
%   [Wrr, Wrf] = generate_weights_random(200,50,180, [0.05 0.02; 0.15 0.1], [0.2; 0.05], 'Seed', 42);
%
% Author: Refactored for clarity & robustness.
% -------------------------------------------------------------------------

p=inputParser; p.FunctionName=mfilename;
addParameter(p,'Seed',[],@(x)isnumeric(x)&&isscalar(x) || isempty(x));
parse(p,varargin{:}); NV=p.Results;
if ~isempty(NV.Seed), rng(NV.Seed,'twister'); end

validateCounts(Ne,Ni,Nx);
validateattributes(Prr, {'double','single'},{'size',[2,2]}, mfilename, 'Prr');
validateattributes(Prx, {'double','single'},{'size',[2,1]}, mfilename, 'Prx');

pee0=Prr(1,1); pei0=Prr(1,2); pie0=Prr(2,1); pii0=Prr(2,2);
pex0=Prx(1);   pix0=Prx(2);

Kee=ceil(pee0*Ne); Kei=ceil(pei0*Ne);
Kie=ceil(pie0*Ni); Kii=ceil(pii0*Ni);
Kex=ceil(pex0*Ne); Kix=ceil(pix0*Ni);

Ke=Kie+Kee; Ki=Kei+Kii; Kx=Kex+Kix;
Wrr=zeros(Ke*Ne+Ki*Ni,1,'int32');
Wrf=zeros(Kx*Nx,1,'int32');

for j=1:Ne
    Wrr((1+(j-1)*Ke):(Kee+(j-1)*Ke)) = sort(randsample(1:Ne, Kee, true));
    Wrr((Kee+1+(j-1)*Ke):(Kee+Kie+(j-1)*Ke)) = sort(randsample(Ne+1:Ne+Ni, Kie, true));
end
offset = Ne*Ke;
for j=1:Ni
    Wrr((offset+1+(j-1)*Ki):(offset+Kei+(j-1)*Ki)) = sort(randsample(1:Ne, Kei, true));
    Wrr((offset+Kei+1+(j-1)*Ki):(offset+Kei+Kii+(j-1)*Ki)) = sort(randsample(Ne+1:Ne+Ni, Kii, true));
end

for j=1:Nx
    Wrf((Kx*(j-1)+1):(Kx*(j-1)+Kex)) = sort(randsample(1:Ne, Kex, true));
    Wrf((Kex+1+(j-1)*Kx):(Kex+Kix+(j-1)*Kx)) = sort(randsample(Ne+1:Ne+Ni, Kix, true));
end

end

% ---- local helper ----
function validateCounts(Ne,Ni,Nx)
for v={'Ne','Ni','Nx'}
    name=v{1}; val=eval(name);
    if ~(isnumeric(val)&&isscalar(val)&&val>0&&mod(val,1)==0)
        error('Argument %s must be a positive integer.', name);
    end
end
end
