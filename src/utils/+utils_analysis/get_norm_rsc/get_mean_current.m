function [inputXr, inputEr, inputIr, inputTr] = get_mean_current(rate, rateX, weight_file, param)

load(weight_file, 'presynIDX','presynIDE','presynIDI');

N = size(rate,1); 
M = size(rate,2);

inputX = zeros(N,2); 
inputE = zeros(N,1); 
inputI = zeros(N,1); 

inputXr = zeros(N,M,2); 
inputEr = zeros(N,M); 
inputIr = zeros(N,M); 

for ii = 1:N
    LocY= mod(presynIDX{ii}-1,100)+1;
    inputX(ii,1) = nnz(LocY<=50);  % area 1
    inputX(ii,2) = nnz(LocY>50);  % area 2
    inputXr(ii,:,1) = sum(rateX(presynIDX{ii}(LocY<=50),1:M),1);  % area 1
    inputXr(ii,:,2) = sum(rateX(presynIDX{ii}(LocY>50),1:M),1);  % area 2
    inputE(ii) = nnz(presynIDE{ii}); 
    inputEr(ii,:) = sum(rate(presynIDE{ii},1:M),1); 
    inputI(ii) = nnz(presynIDI{ii}); 
    inputIr(ii,:) = sum(rate(presynIDI{ii}+param.Ne,1:M),1); 
end
inputXr(1:param.Ne,:,:) = inputXr(1:param.Ne,:,:)*param(1).Jx(1)*1e-3; 
inputXr(1+param.Ne:N,:,:) = inputXr(1+param.Ne:N,:,:)*param(1).Jx(2)*1e-3; 
inputEr(1:param.Ne,:) = inputEr(1:param.Ne,:)*param(1).Jr(1,1)*1e-3; 
inputEr(1+param.Ne:N,:) = inputEr(1+param.Ne:N,:)*param(1).Jr(2,1)*1e-3; 
inputIr(1:param.Ne,:) = inputIr(1:param.Ne,:)*param(1).Jr(1,2)*1e-3; 
inputIr(1+param.Ne:N,:) = inputIr(1+param.Ne:N,:)*param(1).Jr(2,2)*1e-3; 

inputTr = inputEr+inputIr+squeeze(sum(inputXr,3));
