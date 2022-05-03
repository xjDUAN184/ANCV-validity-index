clc;clear;
data = cell2mat(struct2cell(load('.\dataset\Frame.mat '))); 
label = cell2mat(struct2cell(load('.\dataset\FrameGroup.mat')));
K = 10;   SNN_thr = 3;optimal = length(unique(label));
%% 数据标准化
data_max = max(data); 
data_min = min(data);
bre = [];lk = 1;
for j=1:size(data,2) 
    if data_max(j) - data_min(j) <= 0.0001
        bre(lk) = j;
        lk = lk + 1;
        continue;
    else
        data(:,j) = (data(:,j)-data_min(j))/(data_max(j)-data_min(j));
    end
end
data(:,bre) = [];
%%
[N,dim] = size(data);
K_max = ceil(sqrt(N));
k_value = 2;
score = zeros(K_max,2);
while k_value <= K_max
     %% 产生候选结果,这里添加您的聚类算法
    output = kmeans(data,k_value);
    %% 对ANCV指标测试
    [Separation,Compactness,Sep_clu,~] = validity_index(data,output,K,SNN_thr);
    score(k_value,1) = Separation;
    score(k_value,2) = Compactness;
    k_value = k_value + 1;
end
Sep = score(:,1);Com = score(:,2);
%%
index = zeros(K_max,1);
for i = 2:K_max
    index(i) = Sep(i) - Com(i);
end
index(1) = -99999;
[~,res1] = max(index);
if res1 == optimal
    str = ["结果正确，类别数为:",num2str(res1)];
    disp(str);
else
    str = ["结果错误，正确类别数为:",num2str(OPTIMAL(u)),"识别出的类别为:",num2str(res1)];
    disp(str);
end
