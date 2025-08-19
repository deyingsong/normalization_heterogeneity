function [IXr, IEr, IIr, normX, normE, normI] = atttRec_calc_current_norm(weight_name,stats_name,J,optSave)

load(weight_name,'Wrr','Wrf',...
    'seqindE', 'seqindI', 'seqindX');
temp1=load(stats_name);

if isfield(temp1, 'Stats1')
    V1fr = temp1.Stats1.V1fr;
    MTfr = temp1.Stats1.MTfr;
else
    V1fr = temp1.V1fr;
    MTfr = temp1.MTfr;
end

N=2.5e4;Nx=5e3;Ne=2e4;Ni=5e3;

inputXr = zeros(N,3); 
inputEr = zeros(N,3); 
inputIr = zeros(N,3); 

for ii=1:Nx
    ind_temp = Wrf(seqindX(ii)+1:seqindX(ii+1));
    inputXr(ind_temp,:) = inputXr(ind_temp,:)+V1fr(ii,:);
end

for ii=1:Ne
    ind_temp = Wrr(seqindE(ii)+1:seqindE(ii+1));
    inputEr(ind_temp,:) = inputEr(ind_temp,:)+MTfr(ii,:);
end

for ii=1:Ni
    ind_temp = Wrr(seqindI(ii)+1+seqindE(end):seqindI(ii+1)+seqindE(end));
    inputIr(ind_temp,:) = inputIr(ind_temp,:)+MTfr(ii+Ne,:);
end

Jex=J.Jex/sqrt(N);
Jix=J.Jix/sqrt(N);
Jee=J.Jee/sqrt(N);
Jei=J.Jei/sqrt(N);
Jie=J.Jie/sqrt(N);
Jii=J.Jii/sqrt(N);

inputXr(1:Ne,:) = inputXr(1:Ne,:)*Jex*1e-3;
inputXr(Ne+1:N,:) = inputXr(Ne+1:N,:)*Jix*1e-3;
inputEr(1:Ne,:) = inputEr(1:Ne,:)*Jee*1e-3;
inputEr(Ne+1:N,:) = inputEr(Ne+1:N,:)*Jie*1e-3;
inputIr(1:Ne,:) = inputIr(1:Ne,:)*Jei*1e-3;
inputIr(Ne+1:N,:) = inputIr(Ne+1:N,:)*Jii*1e-3;

IXr = inputXr;
IEr = inputEr;
IIr = inputIr;

normX = (IXr(1:Ne,1)+IXr(1:Ne,2))./IXr(1:Ne,3);
normE = (IEr(1:Ne,1)+IEr(1:Ne,2))./IEr(1:Ne,3);
normI = (IIr(1:Ne,1)+IIr(1:Ne,2))./IIr(1:Ne,3);


if optSave
    StatsCurrent.IXr = IXr;
    StatsCurrent.IEr = IEr;
    StatsCurrent.IIr = IIr;
    StatsCurrent.normX = normX;
    StatsCurrent.normE = normE;
    StatsCurrent.normI = normI;
    save_name = strrep(stats_name, 'Stats1', 'StatsCurrent');
    save(save_name, 'StatsCurrent');
end