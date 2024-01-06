function result_K = ANCV_index_py(cluster_name, K, SNN_thr, data_path, data_group_path)
    
    path_len = length(data_path);
    %% load data
    if strcmp(data_path(path_len-2:path_len), 'txt') && isempty(data_group_path)
        % If only the path of the txt file is passed in, the txt file contains data and label 
        % the first column to the second to last column are data, and the last column is label
        data_all = importdata(data_path);
        column = size(data_all, 2);
        data = data_all(:, column-1);
        label = data_all(:, column);
    elseif ~isempty(data_group_path) && strcmp(data_group_path(length(data_group_path)-2:length(data_group_path)), 'mat')
        % If you are passing in .mat file, you need to pass in two files: data_path and data_group_path
        % data_path is data, data_group_path is label
        data = cell2mat(struct2cell(load(data_path))); 
        label = cell2mat(struct2cell(load(data_group_path)));
    end
        
    optimal = length(unique(label));
    %% Data normalization
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
    [N, dim] = size(data);
    K_max = ceil(sqrt(N));
    k_value = 2;
    score = zeros(K_max,2);
    while k_value <= K_max
         %% generate candidate results, add your clustering algorithm here
         if strcmp(cluster_name, 'CTCEHC')
            output = CTCEHC(data, k_value); 
         elseif strcmp(cluster_name, 'NTHC')
            output = NTHC_clustering(data,k_value);
         elseif strcmp(cluster_name, 'kmeans')
            output = kmeans(data,k_value);
         end
        %% compute with ANCV index
        [Separation,Compactness,Sep_clu,~] = validity_index(data,output,K,SNN_thr);
        score(k_value,1) = Separation;
        score(k_value,2) = Compactness;
        k_value = k_value + 1;
    end
    Sep = score(:,1);Com = score(:,2);
    %% Count the ANCV index values from 2 to root N, get the optimal number of clusters.
    index = zeros(K_max,1);
    for i = 2:K_max
        index(i) = Sep(i) - Com(i);
    end
    index(1) = -99999;
    [~,res1] = max(index);
    if res1 == optimal
        str = ["The result is correct, the number of categories is:",num2str(res1)];
        disp(str);
    else
        str = ["The result is wrong, the correct number of categories is:",num2str(optimal),"The categories identified are:",num2str(res1)];
        disp(str);
    end
    result_K = num2str(res1);
end