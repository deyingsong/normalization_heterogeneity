function spiking_sim_Fisher_information_match_indegree(varargin)

% mex EIF1DRFfastslowSyn.c

% run U & A, with same s1;
%  for stim_type='OriMap_gabor_Tseg'
% save spike counts only
% weights in int32
% rand number generator is 'combRecursive'
% rng(Wseed1,'combRecursive')
% connectivity depends on tuning similarity

% RF2D3layer(option, ParamChange)

% param is a struc w/ fields: Ne, Ni, Nx, Jx, Jr, Kx, Kr,
%       gl, Cm, vlb, vth, DeltaT, vT, vl, vre, tref, tausyn, V0, T, dt,
%       maxns, Irecord, Psyn
%   Jx=[Jex; Jix]; Jr=[Jee, Jei; Jie, Jii];
%   Kx=[Kex; Kix]; Kr=[Kee, Kei; Kie, Kii]; % out-degrees are fixed
%   taursyn: syn rise time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   taudsyn: syn decay time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   Psyn(i,j): percentage of synapse j for (X, E, I (i=1,2,3))
%   sigmaRR=[sigmaee, sigmaei; sigmaie, sigmaii];
%   sigmaRX=[sigmaeX; sigmaiX];

%   Wrr is a vector of connections among the recurrent layer, containing postsynaptic cell indices,
%       sorted by the index of the presynaptic cell. The block of postsynaptic cell indices for each presynaptic
%       cell is sorted as excitatory followed by inhibitory cells. I use fixed number of projections Kab to each population.
%       For example, Wrr[j*(Kee+Kie)] to Wrr[j*{Kee+Kie)+Kee-1] are connections from j to E pop and
%       Wrr[j*(Kee+Kie)+Kee] to Wrr[(j+1)*{Kee+Kie)-1] are connections from j to I pop.
%   Wrf is a vector of connections from the feedforward layer to the recurrent layer, sorted by the index of the presynaptic cell.
%       The block of postsynaptic cell indices for each presynaptic cell is sorted as excitatory followed by inhibitory cells.

%   conversion of neuron ID (exc) to (x,y) coordinate in [1, Ne1]x[1, Ne1]:
        % exc. ID [1, Ne], x=ceil(I/Ne1); y=(mod((I-1),Ne1)+1); ID=(x-1)*Ne1+y
        % inh. ID [Ne+1, Ne+Ni], x=ceil((I-Ne)/Ni1); y=(mod((I-Ne-1),Ni1)+1); ID=(x-1)*Ni1+y+Ne;

% sx: spike trains from Layer0
%     sx(1,:) contains spike times.
%     sx(2,:) contains indices of neurons that spike
% s1: spike trains from Layer1
% s2: spike trains from Layer2
%
% data save in filename
% options is a struct w/ fields:
%   'save','CompCorr','plotPopR','fixW','savecurrent','loadfr1','Layer1only'. Default values are 0.
% ParamChange is a cell of 2 columns,
%    the 1st column is variable names and
%    the 2nd column is the values.

% if options.save or options.savecurrent is 1, ParamChange needs to have field 'filename'.
% if options.CompCorr is 1, ParamChange needs to have field 'Nc',e.g. Nc=[500 500];
%      # of neurons to sample from Layer2 & Layer3.
% if options.fixW is 1, ParamChange needs to have field 'Wseed1' & 'Wseed2'.

nVarargs = length(varargin);
switch nVarargs
    case 1
        option = varargin{1};
    case 2
        option = varargin{1};
        ParamChange = varargin{2};
end

if ~isfield(option, 'save') option.save=0; end
% TF = isfield(S,field) returns 1 if field is the name of a field of the 
% structure array S. Otherwise, it returns 0.
if ~isfield(option, 'CompCorr') option.CompCorr=0; end
if ~isfield(option, 'loadfr1') option.loadfr1=1; end
if ~isfield(option, 'fixW') option.fixW=0; end
if ~isfield(option, 'useWfile') option.useWfile=0; end
if ~isfield(option, 'plotPopR') option.plotPopR=0; end
if ~isfield(option, 'savecurrent') option.savecurrent=0; end
if ~isfield(option, 'Layer1only') option.Layer1only=0; end
if ~isfield(option, 'saveRm') option.saveRm=0; end
if ~isfield(option, 'saveSx') option.saveSx=0; end
if ~isfield(option, 'saveS2') option.saveS2=1; end
if ~isfield(option, 'saveParam') option.saveParam=0; end
if ~isfield(option, 'saveW') option.saveW=0; end

if option.save==1
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data')
    end
end
if option.CompCorr==1
    if ~ismember('Nc',ParamChange(:,1))
        error('No Nc (1x2): # of neurons to sample to compute correlations')
    end
end
if option.loadfr1==1
    if ~ismember('fr_fname1',ParamChange(:,1))
        error('No fr_fname')
    end
end
if option.useWfile==1
    if ~ismember('W_fname',ParamChange(:,1))
        error('No W_fname')
    end
end
if option.saveW==1
    if ~ismember('W_fname',ParamChange(:,1))
        error('No W_fname')
    end
end
if option.savecurrent==1
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data')
    end
end

%% define parameters
dim ='2D';

Nx1 = 100;
Ne1 = 200;
Ni1 = 100;

param.Ne = Ne1*Ne1/2;
param.Ni = Ni1*Ni1/2;
param.Nx = Nx1*Nx1/2;
param.N = param.Ne + param.Ni;



% stim_type= 'LocalCorr'; correlation width 'sigmac'
% stim_type='spatialInput';  Gaussian inputs centered at 'center' with width 'sigmac' and mean rate 'rX'
%   center=[.5 .5];
%   sigmac=.15;  % size 1xNstim
% stim_type= 'GlobalCorr'; % cell of size 1xNstim

T=20000; % Total sim time (in msec)

% static currents to MT
inE=0;
inI=0;

% Connection widths
param.sigmaRX = .1*ones(2,1);
param.sigmaRR = .2*ones(2,2);


% number of neurons to record synaptic inputs and voltages from ???
nrecordE0=100;
nrecordI0=100;

% Synaptic time constants
param.taudsyn=[5 100; 5, 100; 8, 100]; % rows: X, E, I, column for different syn types
param.taursyn=[1 2; 1, 2; 1, 2]; % rows: X, E, I
param.Psyn=[.2 .8; 1, 0; 1, 0]; % percentage of diff syn currents


% Connection probabilities (kind of, see use below)
param.Prr=[.01, .04; .03, .04];
param.Prx=[ .05; .05];

% Connection strengths (scaled by sqrt(N) later)
param.Jr=[80 -240; 40, -300];

param.Iapp=[inE;inI];

dt=.05;  % bin size  % 0.01
Tburn=10;   % Burn-in period

% change parameters
if nVarargs==2
    for i=1:size(ParamChange,1)
        eval([ParamChange{i,1} '= ParamChange{i,2};']);
        % eval Execute string with MATLAB expression. eval(s), 
        % where s is a string, causes MATLAB to execute the string 
        % as an expression or statement.
    end
end

if option.loadfr1
    p_stim.fr_fname1=fr_fname1;
end

if size(param.Iapp,2) ~= size(param.Jx,2)
    error('size(param.Iapp,2) should equal size(param.Jx,2), which is the number of parameter sets ')
else
    Np=size(param.Iapp,2);
end
for pid=1:Np
    fprintf('\ninE=%.2f, inI=%.2f\n',param.Iapp(1,pid),param.Iapp(2,pid))
end

clear ParamChange varargin;

%% initialization
maxrate=.05; 
param.maxns=param.N*T*maxrate;
fprintf('\nmaximum average rates to record: %d\n (Hz)',round(maxrate*1e3))

param.dt=dt;
param.T=T;
% EIF neuron paramters
param.gl=[1/15 1/10];  % E, I
param.Cm=[1 1];
param.vlb=[-100 -100];
param.vth=[-10 -10];
param.DeltaT=[2 .5];
param.vT=[-50 -50]; %mV
param.vre=[-65 -65];
param.tref=[1.5 .5];
V0min=param.vre(1);
V0max=param.vT(1);
param.vl=param.Iapp(:,1)'.*[15, 10]-60;
param.V0=(V0max-V0min).*rand(param.N,1)+V0min;
param.Kr=ceil(param.Prr.*[param.Ne, param.Ne; param.Ni,param.Ni]);
param.Kx=ceil(param.Prx.*[param.Ne; param.Ni]);

param.Irecord=[randi(param.Ne,1,nrecordE0), (randi(param.Ni,1,nrecordI0)+param.Ne)];% neuron indice to record synaptic currents and Vm
param.Jr=param.Jr/sqrt(param.N);
param.Jx=param.Jx/sqrt(param.N);

% % Effective connection weights
% q=param.Ne/param.N;
% wrx=(param.Jx(:,1)).*param.Prx*param.Nx/param.N;
% wrr=(param.Jr).*param.Prr.*[q, 1-q; q, 1-q];

%param.Imean=p_stim.rX*wrx*param.N + param.Iapp;
%fprintf('\nfiring rates for large N: %.2f %.2f\n',-(wrr*param.N)\(p_stim.rX*wrx*param.N + param.Iapp)*1e3)

scurr = rng;


%% Simulation

sigma_n=3.5;
tau_n=40;
NI=450;
load('V1filterRecSig0d2Lam0d6.mat');
for i=1:50
    for j=1:50
        for k=1:225
            F1((i-1)*50+j,k)=F((i-1)*100+j,k);
            F2((i-1)*50+j,k)=F((i-1)*100+j+50,k+225);
        end
    end
end



% Simulate Network
if((param.N)<=200000)
    disp('simulation starts')

%     save(filename,'T')
    % Random initial membrane potentials
    if option.loadfr1
        load(fr_fname1, 'fr');
        rng(scurr);
        save(filename,'scurr');
        s1 = attt_genXspk_noise_OnOffCon(fr,T,F1,F2,NI,sigma_n,tau_n,param.ConPerc);
        save(filename,'s1','-append');
        if option.save
            save(filename,'T','-append');
        end
    end
    

    if option.saveParam
        save(filename,'param','p_stim','-append')
    end
    if option.fixW
        if option.useWfile==1
            load(W_fname,'Wrr','Wrf',...
                'outdegreeEE','outdegreeIE','outdegreeEI','outdegreeII','outdegreeEX','outdegreeIX','seqindX','seqindE','seqindI');
            Wrr2 = Wrr; Wrf2 = Wrf;
            param.outdegreeEE = int32(outdegreeEE);
            param.outdegreeEI = int32(outdegreeEI);
            param.outdegreeIE = int32(outdegreeIE);
            param.outdegreeII = int32(outdegreeII);
            param.outdegreeEX = int32(outdegreeEX);
            param.outdegreeIX = int32(outdegreeIX);
            param.seqindX = int32(seqindX);
            param.seqindE = int32(seqindE);
            param.seqindI = int32(seqindI);
            param.Kemax = max(outdegreeEE+outdegreeIE);
            param.Kimax = max(outdegreeEI+outdegreeII);
            param.Kxmax = max(outdegreeEX+outdegreeIX);
            clear Wrr Wrf;
            fprintf('load weight from %s\n',W_fname)
        end
    end
    
    rng(scurr);
    Jx=param.Jx;
    
    tic
    param.vl=param.Iapp(:,1)'.*[15, 10]-60;
    param.V0=(V0max-V0min).*rand(param.N,1)+V0min;
    param.Jx=Jx(:,1);
    [s2,~,~]=EIF1DRFfastslowSynAtttSpatRecInDegree(s1, Wrf2,Wrr2,param);
    elapsetime=toc;

    End=find(s2(2,:)==0,1)-1;
    
    if option.saveRm
        re2=hist(s2(1,s2(2,:)<=param.Ne&s2(2,:)>0),1:T)/param.Ne*1e3;
        ri2=hist(s2(1,s2(2,:)>param.Ne),1:T)/param.Ni*1e3;
        save(filename,'re2','ri2','-append')
    end
    nuSim(1)=1000*nnz(s2(1,1:End)>Tburn & s2(2,1:End)<=param.Ne)/(param.Ne*(T-Tburn));  % Hz
    nuSim(2)=1000*nnz(s2(1,1:End)>Tburn & s2(2,1:End)>param.Ne)/(param.Ni*(T-Tburn));
    fprintf('\naverage rates \n E2: %.2f, I2: %.2f,\n elapsetime=%.2f sec\n',...
        nuSim(1),nuSim(2),elapsetime)

    if option.save
        if option.saveS2
            s2=s2(:,s2(2,:)~=0);
            MTE2 = spktime2count(s2,1:20000,100,200,1);
            timeinds = reshape(((3:40)-1)*5+(1:2)',[],1);
            MTE2_process = MTE2(:,timeinds);
            MTE2_process = MTE2_process(:,1:2:75)+MTE2_process(:,2:2:76);
            MTE2_process = int8(MTE2_process');
            save(filename,'s2','MTE2','MTE2_process','nuSim','-append')
            
        end
    end
    param.Jx=Jx;

    clear Wrr2 Wrf2;
else
    error('N too large') % Your computer probably can't handle this
end
disp('simulation ends')



