clear all

% data = load('Aggregation.txt');
data = cell2mat(struct2cell(load('D:\MATLABprogram\新建文件夹\CciMST-master\data_sets\Glass.mat')));
display(size(data, 1));

MAKER_SIZE = 5;
N = size(data, 1);
clusters = 10;
minN = sqrt(N);

tic
[labels, mst] = CTCEHC(data, clusters, minN, 1, 0);
toc


% cmap = colormap;
% 
% hold on
% for i = 1:clusters
%         ic = uint8((i*64.) / (clusters*1.));
%         plot(data(labels==i, 1), data(labels==i, 2), 'o', 'MarkerSize', ...
%             4, 'MarkerFaceColor', cmap(ic, :), 'MarkerEdgeColor', cmap(ic, :));
% end
