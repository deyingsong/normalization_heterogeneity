file_prefix='Default_AttFar_ori_Ori1_0_Con1_0d5__';

data_folder = '/ix/chuang/des307/';
Ntr = 38000;
N = 2e4;
X = zeros(Ntr,N);

data = load([data_folder file_prefix 'cond1_1_1.mat']);
X(1:Ntr/2,:) = double(data.X1);
data = load([data_folder file_prefix 'cond2_1_1.mat']);
X(Ntr/2+1:Ntr,:) = double(data.X2);

fr = mean(X,1)/0.2;
ind1 = find(fr>2);
save([data_folder file_prefix 'ind1.mat'], 'ind1');