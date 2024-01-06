%%% 画图
data1 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\D4.mat')));  % DS1
data2 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\D2.mat')));  % DS2
data3 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\ls.mat')));  % DS3
data4 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\cth.mat')));  % DS4
data5 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\banana.mat')));  % DS5
data6 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\D12.mat')));  % DS6
data7 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\2d-4c-n4.mat')));  % DS7
data8 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\smileface.mat')));  % DS8
data9 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\donut3.mat')));  % DS9
data10 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\chainlink.mat')));  % DS10
data11 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\lsun.mat')));  % DS11
data12 = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\ANCV_index\dataset\Frame.mat')));  % DS12

data = data8;
%% 标准化
% data_max = max(data); 
% data_min = min(data);
% bre = [];lk = 1;
% for j=1:size(data,2) 
%     if data_max(j) - data_min(j) <= 0.000001
%         bre(lk) = j;
%         lk = lk + 1;
%         continue;
%     else
%         data(:,j) = (data(:,j)-data_min(j))/(data_max(j)-data_min(j));
%     end
% end
% data(:,bre) = [];
%% 聚类
% output = CTCEHC(data, 2);
output = NTHC_clustering(data, 4);

% data_max = max(data); 
% data_min = min(data);
% bre = [];lk = 1;
% for j=1:size(data,2) 
%     if data_max(j) - data_min(j) <= 0.000001
%         bre(lk) = j;
%         lk = lk + 1;
%         continue;
%     else
%         data(:,j) = (data(:,j)-data_min(j))/(data_max(j)-data_min(j));
%     end
% end
% data(:,bre) = [];

%% 画图
plot(data(:, 1), data(:, 2), '.','color','k','MarkerSize',15);
set(gca,'FontName','Times New Roman','FontSize',13);
hold on;  
output = label
plot(data(output==1,1),data(output==1,2),'r.','MarkerSize',15);
hold on;
plot(data(output==2,1),data(output==2,2),'g.','MarkerSize',15);
plot(data(output==3,1),data(output==3,2),'k.','MarkerSize',15);
plot(data(output==4,1),data(output==4,2),'y.','MarkerSize',15);
plot(data(output==5,1),data(output==5,2),'b.','MarkerSize',15);
plot(data(output==6,1),data(output==6,2),'m.','MarkerSize',15);
% set(gca,'xtick',[-22:5:13],'FontName','Times New Roman','FontSize',18);
% set(gca,'ytick',[-12:4:8],'FontName','Times New Roman','FontSize',18);
% xlim([-22 13]);
% ylim([-12 8]);
set(gca,'FontName','Times New Roman','FontSize',18);
set(gca,'linewidth',2)
% plot(data(output==7,1),data(output==7,2),'c.','MarkerSize',15);
% set(gca,'FontName','Times New Roman','FontSize',13);